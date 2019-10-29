import argparse

import vcf as pyvcf
from typing import Dict, Set

from callers.ab_denovo_caller import ABDenovoCaller
from utils.case_utils import parse_fam_file

HEADER_FILE_NAME = "new_calls_header.vcf"
CALLS_FILE_NAME = "new_calls.tsv"
LIMIT = 10000

class Harness():
    def __init__(self, vcf_file: str, family: Dict, callers: Set, flush = None) -> None:
        super().__init__()
        self.vcf_reader = pyvcf.Reader(filename=vcf_file)
        self.family = family
        self.calls = dict()
        self.callers = callers
        if flush is not None and not flush:
            flush = CALLS_FILE_NAME
        self.calls_file = flush
        if self.calls_file:
            self.open_calls(self.calls_file)

        self.ready = False

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
                        calls[tag]="{}={}".format(call[0], call[1])
                    else:
                        calls[tag] = call[0]
            if not calls:
                continue
            chromosome = record.CHROM
            pos = record.POS
            self.calls[(chromosome, pos)] = calls
            if len(self.calls) > LIMIT:
                self.flush_calls(self.calls_file)
                self.calls.clear()

    def get_calls(self):
        return self.calls

    def write_header(self, file_name = None):
        if file_name is None:
            file_name = HEADER_FILE_NAME
        with open(file_name, "w") as header:
            for caller in self.callers:
                header_line = caller.get_header()
                header.write(header_line + '\n')

    def open_calls(self, file_name):
        tags = sorted({caller.get_tag() for caller in self.callers})
        with open(file_name, "w") as f:
            f.write("# CHROM\tPOS\t{}\n".format('\t'.join(tags)))

    def write_calls(self, file_name = None):
        if file_name is None:
            file_name = CALLS_FILE_NAME
        self.open_calls(file_name)
        self.flush_calls(file_name)

    def flush_calls(self, file_name):
        tags = sorted({caller.get_tag() for caller in self.callers})
        with open(file_name, "w") as f:
            f.write("# CHROM\tPOS\t{}\n".format('\t'.join(tags)))
            for key in self.calls:
                calls = self.calls[key]
                p = [str(k) for k in key]
                line = '\t'.join(p + [calls.get(tag, "") for tag in tags])
                f.write(line + '\n')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Run BGM variant callers")
    parser.add_argument("-i", "--input", "--vcf", dest = "vcf", help="Input VCF file", required=True)
    parser.add_argument("-f", "--family", help="Family (fam) file", required=True)
    args = parser.parse_args()
    print (args)

    vcf_file = args.vcf
    fam_file = args.family

    family = parse_fam_file(fam_file)

    callers = {ABDenovoCaller(recall_genotypes=False)}

    harness = Harness(vcf_file, family, callers)
    harness.write_header()
    harness.run()
    harness.write_calls()

    print("All Done")
