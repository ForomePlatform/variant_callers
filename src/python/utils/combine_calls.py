import argparse
import os
import glob
import sortedcontainers
from typing import Dict

from callers.harness import CALLS_FILE_NAME
from utils.tsv import TSVReader


class CallSet:
    def __init__(self, pos) -> None:
        super().__init__()
        self.pos = pos
        self.samples = []

    def __len__(self) -> int:
        return len(self.samples)



def read_all_calls(folder) -> (Dict, str):
    files = glob.glob(os.path.join(folder, "*.tsv"))
    calls = sortedcontainers.SortedDict()
    headerline = "#"
    for f in files:
        with open(f) as ff:
            headerline = ff.readline()
        tsv_reader = TSVReader(f, "", [])
        for pos in tsv_reader.candidate_calls:
            samples = tsv_reader.candidate_calls[pos]
            if pos in calls:
                call_set = calls[pos]
            else:
                call_set = CallSet(pos)
                calls[pos] = call_set
            for s in samples:
                call_set.samples.append(s.strip())
    return calls, headerline


def write_all_calls(output, headerline, calls):
    with open(output, "w") as tsv:
        tsv.write(headerline)
        for pos in calls:
            samples = ','.join(calls[pos].samples)
            line = "{chromosome}\t{pos}\t{samples}\n".\
                format(chromosome=pos.chromosome, pos=pos.pos, samples=samples)
            tsv.write(line)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Run BGM variant callers")
    parser.add_argument("-f", "--folder",
                        help="Folder, containing individual call files, *.tsv",
                        required=True)
    parser.add_argument("-o", "--output",
            help="Output file with combined calls",
            required=False)

    args = parser.parse_args()
    if not args.output:
        args.output = CALLS_FILE_NAME

    calls, headerline = read_all_calls(args.folder)
    write_all_calls(args.output, headerline, calls)