#  Copyright (c) 2019. Sergey Trifonov, Michael Bouzinier, Shamil Sunyaev
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
from callers.harness import Harness
from utils.case_utils import parse_fam_file


def run (args):
    vcf_file = args.vcf
    fam_file = args.family

    family = parse_fam_file(fam_file)
    if (args.dnlib):
        denovo_caller = BayesDenovoCaller(ABDenovoCaller(), args.results, args.dnlib)
    else:
        denovo_caller = ABDenovoCaller()

    callers = {
        denovo_caller,
        ABCompoundHeterozygousCaller(),
        ABHomozygousRecessiveCaller()
    }

    harness = Harness(vcf_file, family, callers, flush=True)
    harness.write_header()
    t = harness.run()
    n = harness.variant_counter

    print("Processed {:d} variants in {:.2f} seconds; rate = {:3.3f} variant/sec".
          format(n, t, n/t))
    print("Detected {:d} calls in {:d} variants".format(harness.call_counter,
                                                       harness.variant_called))

    harness.write_calls()
    harness.apply_calls("xx.vcf")

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
    args = parser.parse_args()
    print(args)

    run(args)

