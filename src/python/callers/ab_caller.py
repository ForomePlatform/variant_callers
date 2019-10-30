#  Copyright (c) 2019. Sergey Trifonov, Michael Bouzinier, Ignat Leshiner,
#  Shamil Sunyaev
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
        if not self.recall_genotypes:
            self.genotypes = {s.sample: s.gt_type for s in record.samples}
            return self.genotypes

        self.genotypes = dict()
        for sample in record.samples:
            gt = sample.gt_type
            ad = getattr(sample.data, 'AD', None)
            if ad:
                n_ref = ad[0]
                n_alt = sum(ad[1:])
                if n_ref + n_alt > 0:
                    balance = n_alt / (n_ref + n_alt)
                    if balance <= self.AB[0]:
                        gt = 0
                    elif balance <= self.AB[1]:
                        gt = 1
                    else:
                        gt = 2
            self.genotypes[sample.sample] = gt
        return self.genotypes

    def get_af(self, genotypes: Dict):
        unrelated_gts = [
            genotypes[id] for id in self.unrelated_samples
                if genotypes[id] is not None
        ]
        if not unrelated_gts:
            return 0.
        return sum(unrelated_gts) / (2. * len(unrelated_gts))

    def validate(self, genotypes: Dict):
        return self.get_af(genotypes) < self.AF_THRESHOLD

    def affected(self, genotypes: Dict) -> List:
        g = [genotypes[id] for id in self.affected_samples]
        return self.gt_to_int(g)

    def unaffected(self, genotypes: Dict):
        g = [genotypes[id] for id in self.unaffected_samples]
        return self.gt_to_int(g)
