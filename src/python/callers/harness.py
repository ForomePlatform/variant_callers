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
import shutil
import time
import traceback
from collections import OrderedDict

import sortedcontainers
import vcf as pyvcf
from vcf.model import _Record
from typing import Dict, Set, List, Collection

from callers.ab_caller import ABCaller
from callers.abstract_caller import AbstractCaller, VariantContext
from utils.vcf_wrappers import JumpVCFReader

HEADER_FILE_NAME = "new_calls_header.vcf"
CALLS_FILE_NAME = "new_calls.tsv"
LIMIT = 100


def execute(cmd):
    print (cmd) 
    os.system(cmd)


def next_chromosome(chromosome:str) -> str:
    prefix = ""
    if chromosome.startswith('chr'):
        chromosome = chromosome[3:]
        prefix = 'chr'
    try:
        c = int(chromosome)
        if c < 22:
           return prefix + str(c+1)
        return prefix + 'X'
    except:
        chromosome = chromosome.upper()
        if chromosome == 'X':
            return prefix + 'Y'
        elif chromosome == 'M':
            return prefix + '1'
        return None


class Harness():
    def __init__(self, vcf_file: str, family: Dict, callers: Set,
                 flush = None, call_set:List = None, start_pos = None,
                 stop = False) -> None:
        super().__init__()
        self.input_vcf = vcf_file
        if call_set:
            self.vcf_reader = JumpVCFReader(filename=self.input_vcf,
                                            call_set=call_set)
        else:
            self.vcf_reader = pyvcf.Reader(filename=self.input_vcf)
        if start_pos:
            x = start_pos.split(':')
            chromosome = x[0].strip()
            if len(x) > 1:
                pos = int(x[1].strip()) - 1
            else:
                pos = None
            print("Jumping to position: {}: {}".format(chromosome, pos))
            self.vcf_reader.fetch(chromosome, pos)
            if stop:
                self.fetch_next = False
            else:
                self.fetch_next = True
        else:
            self.fetch_next = False
        self.family = family
        self.calls = OrderedDict()
        self.callers = callers
        if flush and not isinstance(flush, str):
            flush = CALLS_FILE_NAME
        self.calls_file = flush
        self.calls_file_open = False
        if self.calls_file:
            self.open_calls()
            self.calls_file_open = True

        self.header_file = None
        self.variant_counter = 0
        self.call_counter = 0
        self.variant_called = 0
        self.use_context = len(callers) > 1
        self.shared_context = None
        self.debug_mode = False

    def update_calls(self, caller:AbstractCaller, all_calls: Dict, new_calls: Dict) -> None:
        if (caller.get_n() > 0):
            for c in new_calls:
                value = new_calls[c]
                if (c in all_calls):
                    value = ','.join([all_calls[c], value])
                all_calls[c] = value
        else:
            all_calls.update(new_calls)


    def init_context(self, samples: Set, record: _Record):
        self.shared_context.reset()
        genotypes = ABCaller.calculate_genotypes(record)
        af = ABCaller.calculate_af(genotypes, samples)
        self.shared_context["genotypes"] = genotypes
        self.shared_context["af"] = af

    def run(self):
        t0 = time.time()
        samples = {s for s in self.vcf_reader.samples}
        for caller in self.callers:
            caller.init(self.family, samples)
        if self.use_context:
            self.shared_context = VariantContext()
            for caller in self.callers:
                caller.set_shared_context(self.shared_context)

        #        for record in self.vcf_reader:
        record = None
        while True:
            prev = record
            try:
                record = next(self.vcf_reader)
            except StopIteration:
                if prev and self.fetch_next:
                    chromosome = next_chromosome(prev.CHROM)
                    if chromosome:
                        print("Jumping to {}".format(chromosome))
                        self.vcf_reader.fetch(chromosome)
                    else:
                        break
                else:
                    break

            self.variant_counter += 1
            if (hasattr(self.vcf_reader, "jump")):
                step = 1000
            else:
                step = 10000
            if (self.variant_counter % step) == 0:
                print("Processed {:d} variants in {:7.2f} sec, detected {:d} calls."
                      " Current: {}:{:d}".
                      format(self.variant_counter, time.time() - t0,
                            self.call_counter,
                            record.CHROM,
                            record.POS
                ))

            try:
                if self.use_context:
                    self.init_context(samples, record)
                else:
                    for caller in self.callers:
                        caller.reset_context()
                calls = dict()
                for caller in self.callers:
                    call = caller.make_call(record)
                    if (call):
                        self.update_calls(caller, calls, call)
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
            except Exception as e:
                print("Error in {}: {}".format(record.CHROM, record.POS))
                print(str(e))
                if self.debug_mode:
                    traceback.print_exc()


        if self.calls_file_open and len(self.calls) > 0:
            self.flush_calls()
            print("Totally: processed {:d} variants, flushed {:d} calls".
                  format(self.variant_counter, self.call_counter))

        return (time.time() - t0)

    def get_calls(self):
        return self.calls

    def write_header(self, file_name = None):
        if file_name is None:
            file_name = HEADER_FILE_NAME
        tags = set()
        with open(file_name, "w") as header:
            for caller in self.callers:
                if len(set(caller.get_all_tags()) - tags) == 0:
                    continue
                tags.update(caller.get_all_tags())
                header_line = caller.get_header()
                header.write(header_line + '\n')
        self.header_file = file_name

    def get_tags(self) -> sortedcontainers.SortedSet:
        tags = sortedcontainers.SortedSet()
        for caller in self.callers:
            tags.update(caller.get_all_tags())
        return tags

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
        if not self.calls_file_open:
            self.open_calls()
        if os.path.exists(self.calls_file):
            try:
                shutil.copyfile(self.calls_file, self.calls_file + ".bak")
            except Exception as e:
                print(e)
        self.flush_calls()

    def flush_calls(self):
        if not self.calls:
            return
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

    def apply_calls(self, output_file, tags = None):
        self.flush_calls()
        if not tags:
            tags = [t for t in self.get_tags()]
        calls_file_final = self.calls_file + ".final"
        if os.path.exists(self.calls_file):
            execute("cp {} {}".format(self.calls_file, calls_file_final))
        else:
            execute("cp {} {}".format(calls_file_final, self.calls_file))

        execute("bgzip -f {}".format(self.calls_file))
        execute("tabix -s1 -b2 -e2 -f {}.gz".format(self.calls_file))
        columns = ','.join(["CHROM","POS"] + tags)
        execute("bcftools annotate -a {}.gz -h {} -c {} -o {} {}".
                  format(self.calls_file, self.header_file, columns, output_file, self.input_vcf))


    @classmethod
    def read_header(cls, header_file):
        metadata = dict()
        with open(header_file) as hdr:
            for line in hdr:
                if not line.startswith("##INFO="):
                    continue
                info = line.strip()[len("##INFO=<"):-1]
                id = None
                number = None
                for x in info.split(','):
                    xx = x.split('=')
                    if xx[0] == "ID":
                        id = xx[1]
                    elif xx[0] == "Number":
                        if xx[1] in ['A', 'R', '.']:
                            number = xx[1]
                        else:
                            number = int(xx[1])
                metadata[id] = number
        return metadata