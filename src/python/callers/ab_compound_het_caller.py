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
import vcf as pyvcf

from .ab_caller import ABCaller


class ABCompoundHeterozygousCaller(ABCaller):
    def __init__(self, recall_genotypes = True):
        super().__init__(recall_genotypes)

    def check_genotypes(self, a: List, u: List) -> Tuple:
        if not all([g > 0 for g in a]):
            return ()
        value = sum(1 << i for i, gt in enumerate(u) if gt > 0)
        if (value):
            return (self.get_my_tag(), str(value))
        return ()

    def get_my_tag(self):
        return super(ABCompoundHeterozygousCaller, self).get_my_tag() + "_CMPD_HET"

    def get_type(self):
        return "Integer"

    def get_description(self):
        return "Compound heterozygous by BGM allele balance caller"
