#  Copyright (c) 2019. Sergey Trifonov, Michael Bouzinier
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

import glob
import os
import sortedcontainers

from .misc import raiseException


def parse_fam_file(fam_file):
    samples = sortedcontainers.SortedDict()

    case_dir, file_name = os.path.split(fam_file)
    case = file_name.split('.')[0]
    map_file = None
    maps = glob.glob(os.path.join(case_dir, "samples*"))
    if len(maps) == 1:
        map_file = maps[0]
    elif len(maps) > 1:
        maps = [m for m in maps if case in m]
        if (len(maps) > 0):
            map_file = maps[0]

    with open (fam_file) as input:
        for line in input:
            if line.startswith('#') or not line:
                continue
            try:
                family, id, father, mother, sex, affected = line.split()
                sample = dict()
                sample['family']    = family
                sample['id']        = id
                sample['father']    = father
                sample['mother']    = mother
                sample['sex']       = int(sex)
                sample['affected']  = (int(affected) == 2)
                samples[id] = sample
            except:
                raiseException('Could not parse fam file line: {}'
                                .format(line.strip()))

    return samples
