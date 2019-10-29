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

import argparse
import os

import vcf as pyvcf
from typing import Dict, Set, List

from callers.ab_denovo_caller import ABDenovoCaller
from utils.case_utils import parse_fam_file

HEADER_FILE_NAME = "new_calls_header.vcf"
CALLS_FILE_NAME = "new_calls.tsv"
LIMIT = 10000


def execute(cmd):
    print (cmd) 
    os.system(cmd)

class Harness():
    def __init__(self, vcf_file: str, family: Dict, callers: Set, flush = None) -> None:
        super().__init__()
        self.input_vcf = vcf_file
        self.vcf_reader = pyvcf.Reader(filename=self.input_vcf)
        self.family = family
        self.calls = dict()
        self.callers = callers
        if flush is not None and not flush:
            flush = CALLS_FILE_NAME
        self.calls_file = flush
        self.calls_file_open = False
        if self.calls_file:
            self.open_calls()
            self.calls_file_open = True

        self.ready = False
        self.header_file = None

    def run(self):
        for record in self.vcf_reader:
            if not self.ready:
                for caller in self.callers:
                    samples = {s.sample for s in record.samples}
                    caller.init(self.family, samples)

            calls = dict()
            for caller in self.callers:
                call = caller.make_call(record)
                if (call):
                    tag = caller.get_tag()
                    if (caller.get_n()):
                        calls[tag]=call[1]
                    else:
                        calls[tag] = "1"
            if not calls:
                continue
            chromosome = record.CHROM
            pos = record.POS
            self.calls[(chromosome, pos)] = calls
            if len(self.calls) > LIMIT:
                self.flush_calls()

    def get_calls(self):
        return self.calls

    def write_header(self, file_name = None):
        if file_name is None:
            file_name = HEADER_FILE_NAME
        with open(file_name, "w") as header:
            for caller in self.callers:
                header_line = caller.get_header()
                header.write(header_line + '\n')
        self.header_file = file_name

    def get_tags(self) -> List:
        return sorted({caller.get_tag() for caller in self.callers})

    def open_calls(self):
        tags = self.get_tags()
        with open(self.calls_file, "w") as f:
            f.write("# CHROM\tPOS\t{}\n".format('\t'.join(tags)))
        self.calls_file_open = True

    def write_calls(self, file_name = None):
        if file_name:
            self.calls_file = file_name
        elif not self.calls_file:
            self.calls_file = CALLS_FILE_NAME
        self.open_calls()
        self.flush_calls()

    def flush_calls(self):
        if not self.calls_file_open:
            self.open_calls()
        tags = self.get_tags()
        with open(self.calls_file, "a") as f:
            for key in self.calls:
                calls = self.calls[key]
                p = [str(k) for k in key]
                line = '\t'.join(p + [calls.get(tag, "") for tag in tags])
                f.write(line + '\n')
        self.calls.clear()

    def apply_calls(self, output_file):
        self.flush_calls()
        tags = self.get_tags()
        
        execute("bgzip -f {}".format(self.calls_file))
        execute("tabix -s1 -b2 -e2 -f {}.gz".format(self.calls_file))
        columns = ','.join(["CHROM","POS"] + tags)
        execute("bcftools annotate -a {}.gz -h {} -c {} -o {} {}".
                  format(self.calls_file, self.header_file, columns, output_file, self.input_vcf))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Run BGM variant callers")
    parser.add_argument("-i", "--input", "--vcf", dest = "vcf", help="Input VCF file", required=True)
    parser.add_argument("-f", "--family", help="Family (fam) file", required=True)
    args = parser.parse_args()
    print (args)

    vcf_file = args.vcf
    fam_file = args.family

    family = parse_fam_file(fam_file)

    callers = {ABDenovoCaller(recall_genotypes=True)}

    harness = Harness(vcf_file, family, callers)
    harness.write_header()
    harness.run()
    harness.write_calls()
    harness.apply_calls("xx.vcf")

    print("All Done")
