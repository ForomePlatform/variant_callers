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

from typing import Dict, Set, Tuple, List
import vcf as pyvcf

from .ab_caller import ABCaller


class ABHomozygousRecessiveCaller(ABCaller):
    def __init__(self, recall_genotypes = True):
        super().__init__(recall_genotypes)

    def check_genotypes(self, a: List, u: List) -> Tuple:
        if (all([g == 2 for g in a]) and all([g < 2 for g in u])):
            return (self.get_my_tag(), None)
        return ()

    def get_my_tag(self):
        return super(ABHomozygousRecessiveCaller, self).get_my_tag() + "_HOM_REC"

    def get_type(self):
        return "Flag"

    def get_description(self):
        return "Homozygous recessive by BGM allele balance caller"
