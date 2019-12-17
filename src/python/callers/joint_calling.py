#  Copyright (c) 2019. Partners HealthCare, Harvard Medical School’s
#  Department of Biomedical Informatics
#
#  Developed by Michael Bouzinier, based on contributions by:
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
import vcf as pyvcf

from callers.ab_denovo_caller import SpABDenovoCaller
from callers.harness import Harness
from callers.joint_denovo_caller import JointDenovoCaller
from utils.case_utils import parse_all_fam_files, get_trios_for_family
from utils.tsv import create_tsv_reader


def run(args):
    vcf_file = args.vcf
    # families = parse_all_fam_files(args.families)
    # vcf_reader = pyvcf.Reader(filename=vcf_file)
    # callers = {
    #     SpABDenovoCaller(families[name]) for name in families
    #         if (all([s in vcf_reader.samples for s in families[name]])
    #                and len(get_trios_for_family(families[name])) > 0)
    # }

    families = None
    call_set = None
    calls_file = args.f1
    all_families = parse_all_fam_files(args.families)
    if args.filter:
        if '-' in args.filter:
            x = args.filter.split('-')
            families = [f for f in all_families if x[0] <= f <= x[1]]
        else:
            families = args.filter.split(',')

    if calls_file:
        tsv_reader = create_tsv_reader(families, all_families, calls_file,
                                       "{sample}:PASSED")
        call_set = tsv_reader.call_list()

    callers = {JointDenovoCaller(f_metadata=args.families, vcf_file=vcf_file,
                                 path_to_bams=args.bams, path_to_library=args.dnlib,
                                 bayesian=True, first_stage_calls=calls_file,
                                 families_subset=families)}

    if args.output:
        flush = args.output
    else:
        flush = True

    harness = Harness(vcf_file, family=None, callers=callers, flush=flush,
                      call_set=call_set, start_pos=args.start)
    harness.write_header()
    t = harness.run()
    n = harness.variant_counter

    print("Processed {:d} variants in {:.2f} seconds; rate = {:3.3f} variant/sec".
          format(n, t, n/t))
    print("Detected {:d} calls in {:d} variants".format(harness.call_counter,
                                                       harness.variant_called))

    harness.write_calls()
    if args.apply:
        harness.apply_calls("xx.vcf")

    print("All Done")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Run BGM variant callers")
    parser.add_argument("-i", "--input", "--vcf", dest="vcf",
            help="Input VCF file, required. Use jointly called VCF "
                             "for better results",
            required=True)
    parser.add_argument("-f", "--families",
                help="Path to collection of Family (fam) files, required",
                required=True)
    parser.add_argument("--bams",
            help="Results directory, only used for running Bayesian callers",
            required=False)
    parser.add_argument("--dnlib",
            help="Path to De-Novo library. If specified, then Bayesian De-Novo "
                 "Caller is used, otherwise Allele Balance Caller",
            required=False)
    parser.add_argument("--1", dest="f1",
            help="Path to a file containing results from teh first stage",
            required=False)
    parser.add_argument("--filter",
            help="Comma separated list of families to include or a range",
            required=False)
    parser.add_argument("--start",
            help="Start position in input VCF File",
            required=False)
    parser.add_argument("--output",
            help="Output file with new calls",
            required=False)
    parser.add_argument("--apply", action="store_true",
            help="", required=False)

    args = parser.parse_args()
    print(args)

    run(args)
