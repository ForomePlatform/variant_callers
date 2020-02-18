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

import argparse

from callers.ab_compound_het_caller import ABCompoundHeterozygousCaller
from callers.ab_denovo_caller import ABDenovoCaller
from callers.ab_homo_rec_caller import ABHomozygousRecessiveCaller
from callers.bayes_denovo_caller import BayesDenovoCaller
from callers.harness import Harness, HEADER_FILE_NAME, CALLS_FILE_NAME
from callers.tag_caller import TagCaller
from utils.case_utils import parse_fam_file


def run (args):
    vcf_file = args.vcf
    fam_file = args.family

    family = parse_fam_file(fam_file)
    b = args.callers and ("de-novo-b" in args.callers)

    need_standard_de_novo = (not args.callers) or ("de-novo" in args.callers) or b

    if (args.dnlib and need_standard_de_novo):
        denovo_caller = BayesDenovoCaller(ABDenovoCaller(),
                                          args.results,
                                          args.dnlib,
                                          include_parent_calls=not b,
                                          base_ref=args.assembly)
    else:
        denovo_caller = ABDenovoCaller()

    if args.callers:
        callers = set()
        if {"de-novo", "de-novo-b"} & set(args.callers):
            callers.add(denovo_caller)
        elif "compound_het" in args.callers:
            callers.add(ABCompoundHeterozygousCaller(),)
        elif "homo-rec" in args.callers:
            callers.add(ABHomozygousRecessiveCaller)
        elif "de-novo2" in args.callers:
            callers.add(BayesDenovoCaller(TagCaller("BGM_BAYES_DE_NOVO"),
                        args.results, args.dnlib, include_parent_calls=False))
    else:
        callers = {
            denovo_caller,
            ABCompoundHeterozygousCaller(),
            ABHomozygousRecessiveCaller()
        }

    if args.output and args.execute:
        flush = args.output
    elif not args.execute:
        flush = False
    else:
        flush = True

    harness = Harness(vcf_file, family, callers, flush=flush,
                      start_pos=args.start, stop = args.stop)
    if args.execute:
        harness.write_header()
        t = harness.run()
        n = harness.variant_counter

        print("Processed {:d} variants in {:.2f} seconds; rate = {:3.3f} variant/sec".
              format(n, t, n/t))
        print("Detected {:d} calls in {:d} variants".format(harness.call_counter,
                                                           harness.variant_called))

        harness.write_calls()
        tags = None
    else:
        harness.calls_file = args.output
        harness.header_file = args.header
        tags = [t for t in Harness.read_header(args.header)]
    if args.apply:
        harness.apply_calls("xx.vcf", tags)

    print("All Done")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Run BGM variant callers")
    parser.add_argument("-i", "--input", "--vcf", dest="vcf",
            help="Input VCF file, required. Use jointly called VCF "
                             "for better results",
            required=True)
    parser.add_argument("-f", "--family", help="Family (fam) file, required",
                        required=True)
    parser.add_argument("--results",
            help="Results directory, only used for running Bayesian callers",
            required=False)
    parser.add_argument("--dnlib",
            help="Path to De-Novo library. If specified, then Bayesian De-Novo "
                 "Caller is used, otherwise Allele Balance Caller",
            required=False)
    parser.add_argument("--callers", nargs="*",
                        help="List of callers if different from default")
    parser.add_argument("--start",
            help="Start position in input VCF File",
            required=False)
    parser.add_argument("--stop", action="store_true",
            help="If start position is given then tells if to stop when reaches "
                 "the end of chromosome",
            required=False)
    parser.add_argument("--output",
            help="Output file with new calls",
            required=False)
    parser.add_argument("--assembly", default="hg19",
            help="Assembly to be used: hg19/hg38",
            required=False)
    parser.add_argument("--header",
            help="File with additional VCF headers",
            required=False)
    parser.add_argument("--apply", action="store_true",
            help="Apply calls after execution is complete", required=False)
    parser.add_argument("--apply_calls", help="Do not run the caller, just "
                        "apply calls from a given file", required=False)

    args = parser.parse_args()

    if args.apply_calls:
        args.apply = True
        args.execute = False
        args.output = args.apply_calls
        if not args.header:
            if args.output == CALLS_FILE_NAME:
                args.header = HEADER_FILE_NAME
            else:
                args.header = args.output + ".hdr"
    else:
        args.execute = True

    print(args)

    run(args)

