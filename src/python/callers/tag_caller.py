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

from .ab_caller import ABCaller
from vcf.model import _Record


class TagCaller(ABCaller):
    def __init__(self, tag:str, recall_genotypes = True):
        super().__init__(recall_genotypes)
        self.tag = tag

    def make_call(self, record: _Record) -> Dict:
        if self.tag in record.INFO:
            return {self.tag: record.INFO[self.tag]}
        return {}

    def check_genotypes(self, a: List, u: List) -> Tuple:
        raise Exception("Should be never called")

    def get_my_tag(self):
        return None

    def get_type(self):
        return None

    def get_description(self):
        return "Returns calls tagged by a specified tag {} in the INFO".\
            format(self.tag)
