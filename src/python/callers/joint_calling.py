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
from utils.case_utils import parse_all_fam_files, get_trios_for_family


def run(args):
    families = parse_all_fam_files(args.families)
    vcf_file = args.vcf
    vcf_reader = pyvcf.Reader(filename=vcf_file)
    callers = {
        SpABDenovoCaller(families[name]) for name in families
            if (all([s in vcf_reader.samples for s in families[name]])
                   and len(get_trios_for_family(families[name])) > 0)
    }
    harness = Harness(vcf_file, family=None, callers=callers, flush=True)
    harness.write_header()
    t = harness.run()
    n = harness.variant_counter

    print("Processed {:d} variants in {:.2f} seconds; rate = {:3.3f} variant/sec".
          format(n, t, n/t))
    print("Detected {:d} calls in {:d} variants".format(harness.call_counter,
                                                       harness.variant_called))

    harness.write_calls()
    #harness.apply_calls("xx.vcf")

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

    args = parser.parse_args()
    print(args)

    run(args)
