#  Copyright (c) 2019. Partners HealthCare, Harvard Medical School’s
#  Department of Biomedical Informatics
#
#  Developed by Michael Bouzinier, based on contributions by:
#  Andrew Bjonnes,
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
import glob
import os
import tarfile

import sortedcontainers
from typing import Dict

from .misc import raiseException


def parse_fam_files_content(content, name):
    samples = sortedcontainers.SortedDict()
    for line in content:
        if isinstance(line, bytes):
            line = line.decode('utf-8')
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
            if len(samples) == 0:
                if not sample['affected']:
                    raiseException("First sample in {} is expected to be proband but is unaffected".
                                   format(name))
                sample['proband'] = True
            else:
                sample['proband'] = False
            samples[id] = sample
        except Exception as e:
            raiseException('Could not parse fam file line: {}. {}'
                            .format(line.strip(), e))

    return samples


def get_trios_for_family(family:Dict) -> Dict:
    trios = dict()
    for sample in family:
        mother = family[sample]['mother']
        father = family[sample]['father']
        if (mother == '0'):
            mother = None
        if (father == '0'):
            father = None
        if mother and father:
            trio = [mother, father, sample]
            trios[sample] = trio
    return trios


def parse_fam_file(fam_file):
    with open (fam_file) as input:
        return parse_fam_files_content(input, fam_file)


def find_all_fam_files(path):
    if os.path.isfile(path):
        if path.endswith("tgz"):
            mode = "r:gz"
        elif path.endswith("tar"):
            mode = "r"
        else:
            raise Exception("Unsupported file format: {}".format(path))
        tar = tarfile.open(path, mode)
        return {member.name: tar.extractfile(member)
                for member in tar.getmembers()
                if member.name.endswith(".fam")}
    elif os.path.isdir(path):
        files = glob.glob(os.path.join(path, "**/*.fam"))
        return {os.path.basename(os.path.dirname(f)): f for f in files}


def parse_all_fam_files(path):
    fam_files = find_all_fam_files(path)
    families = dict()
    for name in fam_files:
        f = fam_files[name]
        if not f:
            print ("ERROR: {}".format(name))
            continue
        if (isinstance(f, str)):
            with open(f) as xx:
                content = xx.readlines()
        else:
            content = f.readlines()
        if (not content):
            print ("ERROR: {}".format(name))
            continue

        key = os.path.basename(name)
        if (key.endswith(".fam")):
            key = key[:-4]
        try:
            family = parse_fam_files_content(content, name)
        except Exception as e:
            print("ERROR: " + str(e))
            continue
        families[key] = family

    return families