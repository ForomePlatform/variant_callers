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

import os
import time

import vcf as pyvcf
from typing import Dict, Set, List

HEADER_FILE_NAME = "new_calls_header.vcf"
CALLS_FILE_NAME = "new_calls.tsv"
LIMIT = 100


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
        if flush and not isinstance(flush, str):
            flush = CALLS_FILE_NAME
        self.calls_file = flush
        self.calls_file_open = False
        if self.calls_file:
            self.open_calls()
            self.calls_file_open = True

        self.ready = False
        self.header_file = None
        self.variant_counter = 0
        self.call_counter = 0
        self.variant_called = 0

    def run(self):
        t0 = time.time()
        for record in self.vcf_reader:
            if not self.ready:
                for caller in self.callers:
                    samples = {s.sample for s in record.samples}
                    caller.init(self.family, samples)
                self.ready = True
            self.variant_counter += 1

            calls = dict()
            for caller in self.callers:
                call = caller.make_call(record)
                calls.update(call)
            if not calls:
                continue
            chromosome = record.CHROM
            pos = record.POS
            self.calls[(chromosome, pos)] = calls
            self.call_counter += len(calls)
            self.variant_called += 1
            if self.calls_file_open and len(self.calls) > LIMIT:
                self.flush_calls()
                print("Processed {:d} variants, flushed {:d} calls".
                      format(self.variant_counter, self.call_counter))
        return (time.time() - t0)

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
        tags = []
        for caller in self.callers:
            tags.extend(caller.get_all_tags())
        return sorted(tags)

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
                calls = {tag: v if v != None else "1" for tag, v in self.calls[key].items()}
                p = [str(k) for k in key]
                line = '\t'.join(p + [calls.get(tag, ".") for tag in tags])
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


