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

from typing import Dict, Set, Tuple, List
from vcf.model import _Record

from callers.ab_caller import ABCaller
from callers.abstract_caller import AbstractCaller
from denovo2.detect.detect2 import DenovoDetector, VariantHandler
from utils.misc import raiseException


class JointDenovoCaller(AbstractCaller):
    def __init__(self, families: Dict, path_to_bams: str,
                 path_to_library: str,
                 pp_threshold: float = 0.7, include_parent_calls: bool = True):
        super().__init__()
        self.families = None
        self.path_to_bams = path_to_bams
        self.path_to_library = path_to_library
        self.pp_threshold = pp_threshold
        self.detector = None
        self.return_parent_calls = include_parent_calls

    def init(self, family: Dict, samples: Set):
        super(JointDenovoCaller, self).init(family, samples)
        self.parent.init(family, samples)
        trio = self.get_trio()
        if (not trio):
            raiseException(
                "For {} family must contain a trio".format(self.get_my_tag()))
        list_of_bam_files = [
            self.path_to_bams.format(sample) for sample in trio
        ]
        self.detector = DenovoDetector(self.path_to_library,
                                       trio_list=list_of_bam_files)

    def make_call(self, record: _Record) -> Dict:
        result = dict()
        parent_call = self.parent.make_call(record)
        if not parent_call:
            return result
        if (self.return_parent_calls and self.parent.get_type()):
            result.update(parent_call)
        chromosome = record.CHROM
        if not chromosome.startswith("chr"):
            chromosome = "chr" + chromosome
        pos = record.POS
        genotypes = self.parent.get_genotypes(record)
        af = self.parent.get_af(genotypes)
        variant = VariantHandler(chromosome, pos, record.REF, record.ALT, af)
        passed = self.detector.detect(variant)
        if (passed):
            pp = variant.getProp("PP")
            if (pp > self.pp_threshold):
                result[self.get_my_tag()] = str(pp)
        return result

    def get_my_tag(self):
        #return super(JointDenovoCaller, self).get_my_tag() + "_BAYES_DE_NOVO"
        return "BGM_BAYES_DE_NOVO"

    def get_type(self):
        return "String"

    def get_description(self):
        return "Probability of de novo by BGM Bayes caller"

    def close(self):
        if (self.detector):
            self.detector.close()
