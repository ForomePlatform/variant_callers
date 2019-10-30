#  Copyright (c) 2019. Sergey Trifonov, Michael Bouzinier, Ignat Leshiner, Shamil Sunyaev
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


class BayesDenovoCaller(AbstractCaller):
    def __init__(self, parent:ABCaller, path_to_bams:str, path_to_library:str):
        super().__init__()
        self.parent = parent
        self.path_to_bams = path_to_bams
        self.path_to_library = path_to_library
        self.detector = None

    def init(self, family: Dict, samples: Set):
        super(BayesDenovoCaller, self).init(family, samples)
        self.parent.init(family, samples)
        trio = self.get_trio()
        if (not trio):
            raiseException("For {} family must contain a trio".format(self.get_my_tag()))
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
            result[self.get_my_tag()] = str(variant.getProp("PP"))
        return result

    def get_my_tag(self):
        return super(BayesDenovoCaller, self).get_my_tag() + "_BAYES_DE_NOVO"

    def get_all_tags(self):
        return [self.get_my_tag(), self.parent.get_my_tag()]

    def get_type(self):
        return "Float"

    def get_description(self):
        return "Probability of de novo by BGM Bayes caller"

    def get_header(self):
        h1 = super(BayesDenovoCaller, self).get_header()
        h2 = self.parent.get_header()
        return '\n'.join([h1, h2])

    def close(self):
        if (self.detector):
            self.detector.close()