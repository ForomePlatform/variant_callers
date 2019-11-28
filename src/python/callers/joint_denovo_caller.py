#  Copyright (c) 2019. Partners HealthCare, Harvard Medical School’s
#  Department of Biomedical Informatics, Sergey Trifonov
#
#  Developed by Sergey Trifonov and Michael Bouzinier, based on contributions by:
#  Anwoy Kumar Mohanty, Andrew Bjonnes,
#  Ignat Leshchiner, Shamil Sunyaev and other members of Division of Genetics,
#  Brigham and Women's Hospital
#
#  Licensed under the Apache License, Version 2.0 (the "License");
#  you may not use this file except in compliance with the License.
#  You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
#  Unless required by applicable law or agreed to in writing, software
#  distributed under the License is distributed on an "AS IS" BASIS,
#  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#  See the License for the specific language governing permissions and
#  limitations under the License.
import os
from typing import Dict, Set, Tuple, List
from vcf.model import _Record

from callers.ab_denovo_caller import ABDenovoCaller
from callers.abstract_caller import AbstractCaller
from denovo2.detect.detect2 import DenovoDetector, VariantHandler
from utils.misc import raiseException
import sortedcontainers
import vcf as pyvcf

from utils.case_utils import parse_all_fam_files, get_trios_for_family, get_bam_patterns

class LocalCaller:
    def __init__(self, proband:str, caller: ABDenovoCaller, detector: DenovoDetector) -> None:
        super().__init__()
        self.proband = proband
        self.parent_caller = caller
        self.detector = detector


class JointDenovoCaller(AbstractCaller):
    def __init__(self, f_metadata:str, vcf_file:str, path_to_bams: str,
                 path_to_library: str,
                 pp_threshold: float = 0.7, bayesian: bool = True):
        super().__init__()
        self.local_callers = dict()
        self.path_to_bams = path_to_bams
        self.path_to_library = path_to_library
        self.pp_threshold = pp_threshold
        self.bayesian = bayesian
        self.calculates_pp = True if self.path_to_bams else False
        if self.calculates_pp:
            self.format = "{sample}:{pp:1.2f}"
        else:
            self.format = "{sample}:{passed}"
        self.check_families(f_metadata, vcf_file)

    def check_families(self, f_metadata:str, vcf_file:str):
        self.local_callers = sortedcontainers.SortedDict()
        families = parse_all_fam_files(f_metadata)
        vcf_reader = pyvcf.Reader(filename=vcf_file)
        patterns = get_bam_patterns()
        bam_pattern = None
        if self.bayesian and self.calculates_pp:
            bam_pattern = os.path.join(self.path_to_bams, patterns[0])
        samples = {s for s in vcf_reader.samples}

        shared_detector = None
        if not self.calculates_pp:
            shared_detector = DenovoDetector(self.path_to_library)

        for name in families:
            family = families[name]
            if not all([s in samples for s in family]):
                continue
            trios = get_trios_for_family(family)
            for proband in trios:
                trio = trios[proband]
                if self.bayesian:
                    if self.calculates_pp:
                        list_of_bam_files = [
                            bam_pattern.format(sample=sample) for sample in trio
                        ]
                        if not all (os.path.exists(bam) for bam in list_of_bam_files):
                            continue
                        detector = DenovoDetector(self.path_to_library,
                                        trio_list=list_of_bam_files)
                    else:
                        detector = shared_detector
                else:
                    detector = None
                ab_caller = ABDenovoCaller()
                ab_caller.set_shared_context(self.variant_context)
                ab_caller.init(family, samples)
                local_caller = LocalCaller(proband, ab_caller, detector)
                self.local_callers[proband] = local_caller
        return

    def init(self, families: Dict, samples: Set):
        return

    def make_call(self, record: _Record) -> Dict:
        result = []

        self.variant_context.reset()
        samples = {s.sample for s in record.samples}
        genotypes = ABDenovoCaller.calculate_genotypes(record)
        af = ABDenovoCaller.calculate_af(genotypes, samples)
        chromosome = record.CHROM
        if not chromosome.startswith("chr"):
            chromosome = "chr" + chromosome
        pos = record.POS
        self.variant_context["genotypes"] = genotypes
        self.variant_context["af"] = af

        for proband in self.local_callers:
            caller = self.local_callers[proband]
            parent_call = caller.parent_caller.make_call(record)
            if not parent_call:
                continue

            value = None
            if self.bayesian:
                variant = VariantHandler(chromosome, pos, record.REF, record.ALT, af)
                passed = caller.detector.detect(variant)
                if (passed):
                    if caller.detector.gives_pp():
                        pp = variant.getProp("PP")
                        if (pp > self.pp_threshold):
                            value = self.format.format(sample=proband, pp=pp)
                    else:
                        value = self.format.format(sample=proband, passed="PASSED")
            else:
                value = self.format.format(sample=proband, pp=1 - af)
            if (value != None):
                result.append(value)

        if result:
            return {self.get_my_tag(): ','.join(result)}
        return dict()

    def get_my_tag(self):
        #return super(JointDenovoCaller, self).get_my_tag() + "_BAYES_DE_NOVO"
        return "BGM_BAYES_DE_NOVO"

    def get_type(self):
        return "String"

    def get_description(self):
        return "Probability of de novo by BGM Bayes caller"

    def close(self):
        for caller in self.local_callers.values():
            if caller.detector:
                caller.detector.close()

