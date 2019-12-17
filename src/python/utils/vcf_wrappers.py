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

from typing import List
from vcf.parser import Reader

from utils.case_utils import parse_all_fam_files
from utils.tsv import create_tsv_reader


class JumpVCFReader(Reader):
    def __init__(self, call_set: List, fsock=None, filename=None,
                 compressed=None, prepend_chr=False, strict_whitespace=False,
                 encoding='ascii'):
        super().__init__(fsock, filename, compressed, prepend_chr,
                         strict_whitespace, encoding)
        self.call_set = call_set
        self.pos_in_call_set = 0

    def jump(self):
        if (self.pos_in_call_set >= len(self.call_set)):
            return None
        call = self.call_set[self.pos_in_call_set]
        self.pos_in_call_set += 1
        return self.fetch(call.chromosome, call.pos-1, call.pos)

    def __next__(self):
        me = self.jump()
        _record = super().__next__()
        return _record


if __name__ == '__main__':
    '''Test me'''
    families = ['udn0013', 'udn0028']
    metadata = parse_all_fam_files('udn_cases.tgz')
    vcf = "udn0233_wes_final.vep.vcf.gz"
    calls_file = "denovo/denovo_1.tsv"

    tsv_reader = create_tsv_reader(families, metadata, calls_file, "{sample}:PASSED")
    samples = tsv_reader.call_list()

    vcf_reader = JumpVCFReader(call_set=samples, filename=vcf)
    for record in vcf_reader:
        print(record.CHROM, ": ", record.POS)


