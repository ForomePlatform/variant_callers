#  Copyright (c) 2019. Partners HealthCare, Harvard Medical School’s
#  Department of Biomedical Informatics, Sergey Trifonov
#
#  Developed by Michael Bouzinier, based on contributions by:
#  Shamil Sunyaev and other members of Division of Genetics,
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
import os

import boto3
import sortedcontainers
import vcf as pyvcf

from utils.case_utils import parse_all_fam_files, get_trios_for_family

BUCKET = "udn-joint-calling"
PREFIX = "cases/wes/"
pattern = "{sample}.recal.realign.dedup.bam"
patterns = [pattern, pattern + ".bai"]


def download_bams_for_trios(metadata:str, vcf_file:str, key:str, secret:str, dest:str):
    s3 = boto3.client('s3', aws_access_key_id=key, aws_secret_access_key=secret)
    families = parse_all_fam_files(metadata)
    vcf_reader = pyvcf.Reader(filename=vcf_file)

    files = sortedcontainers.SortedDict()
    for name in families:
        family = families[name]
        if not all([s in vcf_reader.samples for s in family]):
            continue
        trios = get_trios_for_family(family)
        for trio in trios.values():
            for sample in trio:
                for p in patterns:
                    file_name = p.format(sample=sample)
                    object_name = PREFIX + sample + '/' + file_name
                    files[file_name] = object_name
    print (files)
    print (len(files))

    if not os.path.isdir(dest):
        print("ERROR: destination directory does not exist {}".format(dest))
        return

    for f in files:
        object_name = files[f]
        target = os.path.join(dest, f)
        print("Downloading: {} ==> {}".format(object_name, target))
        try:
            s3.download_file(BUCKET, object_name, target)
        except Exception as e:
            print("ERROR: " + str(e))

    print("All Done.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run BGM variant callers")
    parser.add_argument("-i", "--input", "--vcf", dest="vcf",
            help="Input VCF file, required. Use jointly called VCF "
                             "for better results",
            required=True)
    parser.add_argument("-f", "--families",
                help="Path to collection of Family (fam) files, required",
                required=True)
    parser.add_argument("--key", help="aws_access_key_id, required", required=False)
    parser.add_argument("--secret", help="aws_secret_access_key, required", required=False)
    parser.add_argument("--dest", help="Destination directory, required", required=True)

    args = parser.parse_args()
    print(args)

    if args.key:
        key = args.key
        secret = args.secret
    else:
        key = None
        secret = None
        configuration = "default"
        section = '[' + configuration + ']'
        f = os.path.join(os.path.expanduser('~'), ".aws", "credentials")
        with open(f) as credentials:
            state = 1
            for line in credentials:
                if state == 1:
                    if line.strip() == section:
                        state = 2
                    continue
                if state == 2:
                    text = line.strip()
                    if text.startswith('['):
                        state = 3
                        break
                    if "aws_access_key_id" in text:
                        key = text.split('=')[1].strip()
                    elif "aws_secret_access_key" in text:
                        secret = text.split('=', 1)[1].strip()

    download_bams_for_trios(args.families, args.vcf, key, secret, args.dest)

