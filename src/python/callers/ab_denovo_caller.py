#  Copyright (c) 2019. Sergey Trifonov, Michael Bouzinier, Ignat Leshiner,
#  Shamil Sunyaev
#
#  Licensed under the Apache License, Version 2.0 (the "License");
#  you may not use this file except in compliance with the License.
#  You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#

from typing import Dict, Set, Tuple, List

from .ab_caller import ABCaller


class ABDenovoCaller(ABCaller):
    def __init__(self, recall_genotypes = True):
        super().__init__(recall_genotypes)

    def check_genotypes(self, a: List, u: List) -> Tuple:
        if (all([g > 0 for g in a]) and sum(u) == 0):
            return (self.get_my_tag(), None)
        return ()

    def get_my_tag(self):
        return super(ABDenovoCaller, self).get_my_tag() + "_DE_NOVO"

    def get_type(self):
        return "Flag"

    def get_description(self):
        return "De novo by BGM allele balance caller"
