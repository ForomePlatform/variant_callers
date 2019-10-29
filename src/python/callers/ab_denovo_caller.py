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

from typing import Dict, Set, Tuple
import vcf as pyvcf

from .ab_caller import ABCaller


class ABDenovoCaller(ABCaller):
    def __init__(self, recall_genotypes = True):
        super().__init__(recall_genotypes)

    def make_call(self, record: pyvcf.Reader) -> Tuple:
        genotypes = self.get_genotypes(record)
        if (not self.validate(genotypes)):
            return ()
        a = self.affected(genotypes)
        u = self.unaffected(genotypes)
        if (all([g > 0 for g in a]) and sum(u) == 0):
            return (self.get_tag(), None)
        return ()

    def get_tag(self):
        return super(ABDenovoCaller, self).get_tag() + "_DE_NOVO"

    def get_type(self):
        return "Flag"

    def get_description(self):
        return "De novo by BGM allele balance caller"
