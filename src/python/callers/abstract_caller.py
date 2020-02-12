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
import sys
from abc import ABC, abstractmethod
from typing import Dict, Set, Tuple, List
from vcf.model import _Record

from utils.misc import raiseException
import random

class VariantContext(dict):
    def __init__(self):
        super().__init__()
        self.random = random.randint(0, sys.maxsize)

    def reset(self):
        self.clear()

    def put(self, key, value):
        self[key] = value


class AbstractCaller(ABC):
    def __init__(self):
        self.family = None
        self.samples = None
        self.variant_context = VariantContext()
        self.shared_context = False
        self.unrelated_samples = set()
        return

    def init(self, family: Dict, samples: Set):
        if (family):
            self.family = family
        if (not all([s in samples for s in self.family])):
            raiseException("Samples {} are not in VCF".format(','.join(self.family)))
        self.samples = samples
        self.unrelated_samples = self.samples - self.family.keys()
        self.affected_samples = [id for id in self.family if self.family[id]['affected']]
        self.unaffected_samples = [id for id in self.family if not self.family[id]['affected']]
        return

    @abstractmethod
    def make_call(self, record: _Record) -> Dict:
        return {}

    def get_my_tag(self):
        return "BGM"

    def get_all_tags(self):
        return [self.get_my_tag()]

    def get_type(self):
        return "String"

    def get_n(self):
        type = self.get_type()
        if type == "Flag":
            return 0
        return 1

    def get_trio(self):
        for sample in self.family:
            s = self.family[sample]
            if (s["proband"]):
                return [s["mother"], s["father"], s["id"]]
        return None

    @abstractmethod
    def get_description(self):
        return ""

    def get_header(self):
        pattern = '##INFO=<ID={tag},Number={n},Type={type},Description="{desc}">'
        return pattern.format(tag=self.get_my_tag(), type=self.get_type(),
                              n=self.get_n(), desc=self.get_description())

    @staticmethod
    def gt_to_int(gl: List):
        return [g if g != None else 0 for g in gl]

    def add_to_context(self, key, value):
        self.variant_context.put(key, value)

    def reset_context(self):
        self.variant_context.reset()

    def set_shared_context(self, shared_ctx: VariantContext):
        self.variant_context = shared_ctx
        self.shared_context = True
