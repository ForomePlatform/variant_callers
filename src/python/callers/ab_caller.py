#  Copyright (c) 2019. Partners HealthCare, Harvard Medical School’s
#  Department of Biomedical Informatics
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

from abc import abstractmethod
from typing import Dict, Set, Tuple, List
import vcf as pyvcf
from vcf.model import _Record

from .abstract_caller import AbstractCaller


class ABCaller(AbstractCaller):
    AB = [0.05, 0.85]
    AF_THRESHOLD = 0.1
    def __init__(self, recall_genotypes: bool):
        super().__init__()
        self.recall_genotypes = recall_genotypes
        self.genotypes = None

    def make_call(self, record: _Record) -> Dict:
        self.genotypes = None
        genotypes = self.get_genotypes(record)
        if (not self.validate(genotypes)):
            return {}
        a = self.affected(genotypes)
        u = self.unaffected(genotypes)
        result = self.check_genotypes(a, u)
        if len(result) > 1:
            return {result[0]: result[1]}
        if result:
            return {result[0]: None}
        return {}

    @abstractmethod
    def check_genotypes(self, a: List, u: List) -> Tuple:
        pass

    def get_genotypes(self, record: _Record):
        if self.genotypes:
            return self.genotypes
        if "genotypes" in self.variant_context:
            self.genotypes = self.variant_context["genotypes"]
            return self.genotypes
        if not self.recall_genotypes:
            self.genotypes = {s.sample: s.gt_type for s in record.samples}
            return self.genotypes

        self.genotypes = self.calculate_genotypes(record)
        return self.genotypes

    @classmethod
    def calculate_genotypes(cls, record: _Record) -> Dict:
        genotypes = dict()
        for sample in record.samples:
            gt = sample.gt_type
            ad = getattr(sample.data, 'AD', None)
            if ad and isinstance(ad, list):
                n_ref = ad[0]
                n_alt = sum(ad[1:])
                if n_ref + n_alt > 0:
                    balance = n_alt / (n_ref + n_alt)
                    if balance <= cls.AB[0]:
                        gt = 0
                    elif balance <= cls.AB[1]:
                        gt = 1
                    else:
                        gt = 2
            genotypes[sample.sample] = gt
        return genotypes

    def get_af(self, genotypes: Dict):
        if "af" in self.variant_context:
            return self.variant_context["af"]

        return self.calculate_af(genotypes, self.unrelated_samples)

    @classmethod
    def calculate_af(cls, genotypes: Dict, samples) -> float:
        gts = [genotypes[id] for id in samples
                if genotypes[id] is not None
        ]
        if not gts:
            return 0.
        return sum(gts) / (2. * len(gts))

    def validate(self, genotypes: Dict):
        return self.get_af(genotypes) < self.AF_THRESHOLD

    def affected(self, genotypes: Dict) -> List:
        g = [genotypes[id] for id in self.affected_samples]
        return self.gt_to_int(g)

    def unaffected(self, genotypes: Dict):
        g = [genotypes[id] for id in self.unaffected_samples]
        return self.gt_to_int(g)

    # @classmethod
    # def liftover_38_to_19(cls):
