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
from functools import total_ordering
from typing import Collection, List, Dict

import sortedcollections

from utils.case_utils import get_trios_for_family

@total_ordering
class Pos:
    def __init__(self, chromosome, pos) -> None:
        self.chromosome = chromosome
        self.pos = pos

    def __str__(self) -> str:
        return "{}:{:d}".format(self.chromosome, self.pos)

    def __eq__(self, o: object) -> bool:
        if not isinstance(o, Pos):
            return False
        return self.chromosome == o.chromosome and self.pos == o.pos

    def __lt__(self, o):
        if not isinstance(o, Pos):
            raise TypeError
        chromosomes = []
        for c in [self.chromosome, o.chromosome]:
            if c.startswith('chr'):
                c = c[3:]
            c = c.lower()
            if c == 'm':
                c = 0
            elif c == 'x':
                c = 24
            elif c == 'y':
                c = 25
            else:
                c = int(c)
            chromosomes.append(c)
        if chromosomes[0] == chromosomes[1]:
            return self.pos < o.pos
        return chromosomes[0] < chromosomes[1]

    def __hash__(self) -> int:
        return hash((self.chromosome, self.pos))


class TSVReader:
    def __init__(self, tsv_calls_file:str, format_string:str, samples:Collection) -> None:
        super().__init__()
        self.samples = samples
        patterns = dict()
        patterns.update({s:format_string.format(sample=s) for s in self.samples})
        self.candidate_calls = sortedcollections.SortedDict()
        with open(tsv_calls_file) as calls:
            for call in calls:
                if call.startswith('#'):
                    continue
                samples=set()
                for sample in patterns:
                    if patterns[sample] in call:
                        samples.add(sample)
                if samples:
                    data = call.split()
                    key = Pos(data[0], int(data[1]))
                    self.candidate_calls[key] = samples

    def has_call(self, chromosome, pos):
        return Pos(chromosome, pos) in self.candidate_calls

    def has_sample(self, chromosome, pos, sample):
        call = self.candidate_calls.get(Pos(chromosome, pos))
        if call:
            return sample in call
        return False

    def call_list(self):
        return sorted([call for call in self.candidate_calls])

    def call_set(self):
        return {call for call in self.candidate_calls}


def create_tsv_reader(families:List,
                      metadata:Dict,
                      tsv_calls_file:str,
                      format_string:str) -> TSVReader:
    samples = set()
    if not families:
        families = [f for f in metadata]
    for f in families:
        family = metadata[f]
        trios = get_trios_for_family(family)
        for s in trios:
            samples.add(s)
    return TSVReader(tsv_calls_file, format_string, samples)


