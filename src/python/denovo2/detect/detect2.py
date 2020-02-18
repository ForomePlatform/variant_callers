#  Copyright (c) 2019. Partners HealthCare, Harvard Medical School’s
#  Department of Biomedical Informatics, Sergey Trifonov
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


from .dn_model import DeNovo_Model, DeNovo_MDL_Reader
from .read_pysam import PysamList, AD_LibCollection

#========================================
# Chromosomes name/num
#========================================
def chrom2num(chrom_name):
    if chrom_name.startswith("chr"):
        chrom_name = chrom_name[3:]
    if chrom_name.isdigit():
        chrom_num = int(chrom_name)
        assert 0 <= chrom_num <= 22
    else:
        chrom_num = {"M": 0, "MT": 0, "X": 23, "Y": 24}.get(chrom_name, 25)
    assert chrom_num < 25, "Bad chromosome: " + chrom_name
    return chrom_num

def num2chrom(chrom_num):
    return (str(chrom_num) if 0 < chrom_num <= 22
        else {0: "MT", 23: "X", 24: "Y"}[chrom_num])

#========================================
# Variant works
#========================================
class VariantHandler:
    def __init__(self, chrom_name, pos, ref, alt, af,
            chrom_num = None, ref_base = None):
        if chrom_num is not None:
            assert chrom_name is None
            chrom_name = "chr" + num2chrom(chrom_num)
        self.mChromName = chrom_name
        self.mChromNum = chrom2num(self.mChromName)
        self.mPos = int(pos)
        self.mRef = ref
        self.mAlt = alt
        self.mAF = float(af)
        self.mProps = dict()
        self.mRefBase = "hg19" if ref_base is None else ref_base

    # Attention: set level = 2 for keep only snip variants
    def isGood(self, level = 1):
        if '*' in self.mAlt or ',' in self.mAlt:
            return False
        if level > 1:
            return len(self.mRef) == 1 and len(self.mAlt) == 1
        return True

    def getChromName(self):
        return self.mChromName

    def getChromNum(self):
        return self.mChromNum

    def getPos(self):
        return self.mPos

    def getRef(self):
        return self.mRef

    def getAlt(self):
        return self.mAlt

    def getAF(self):
        return self.mAF

    def getProp(self, key):
        return self.mProps.get(key)

    def setProp(self, key, val):
        self.mProps[key] = val

    def getRefBase(self):
        return self.mRefBase

    @classmethod
    def parse(cls, line):
        params = line.strip().split()
        chrom_name, pos = params[:2]
        ref, alt = params[3:5]
        af = None
        for pair in params[7].split(';'):
            arg, _, val = pair.partition('=')
            if arg == "ExAC_AF_computed":
                af = float(val)
                break
        assert af is not None
        return VariantHandler(chrom_name, pos, ref, alt, af)

    @classmethod
    def loadFile(cls, filename):
        variants = []
        with open(filename, "r") as inp:
            for line in inp:
                var = cls.parse(line)
                if var.isGood():
                    variants.append(var)
        variants.sort(key = lambda var: (var.getChromNum(), var.getPos()))
        return variants

#========================================
class DenovoDetector:
    def __init__(self, unrelated_dir, trio_filename = None, trio_list = None,
            dump_filename = None, mdl_file = None):
        if trio_filename is not None or trio_list is not None:
            self.mTrioSamFiles = PysamList(
                file_with_list_of_filenames = trio_filename,
                list_of_bams = trio_list)
        else:
            self.mTrioSamFiles = None
    # def __init__(self, unrelated_dir, trio_filename=None, trio_list=None,
    #         dump_filename=None, mdl_file=None):
    #     self.mTrioSamFiles = PysamList(
    #         file_with_list_of_filenames = trio_filename,
    #         list_of_bams = trio_list)
        self.mUnrelLib = None
        self.mDumpFName = dump_filename
        if unrelated_dir is not None:
            self.mUnrelLib = AD_LibCollection(unrelated_dir, self.mDumpFName)
        else:
            assert self.mDumpFName is None
            self.mUnrelMdl = DeNovo_MDL_Reader(mdl_file)

    def close(self):
        if self.mUnrelLib is not None:
            self.mUnrelLib.finishUp()
        else:
            self.mUnrelMdl.close()

    def gives_pp(self):
        return self.mTrioSamFiles is not None

    def detect(self,  variant):
        if self.mDumpFName:
            self.mUnrelLib.mineAD(variant)
        if self.mUnrelLib is not None:
            ad_model = DeNovo_Model.createByADLib(variant, self.mUnrelLib)
        else:
            ad_model = self.mUnrelMdl.getPosModel(variant)
        if self.mTrioSamFiles is None:
            variant.setProp("PASSED", not ad_model.isBad())
        else:
            ad_model.evalVariant(variant, self.mTrioSamFiles)
        return variant.getProp("PASSED")

#========================================
def runner(outfilename, initial_filename, unrelated_dir, mdl_file,
        trio_filename, dump_filename):
    detector = DenovoDetector(unrelated_dir, trio_filename,
        dump_filename = dump_filename, mdl_file = mdl_file)

    variants = VariantHandler.loadFile(initial_filename)

    for count, variant in enumerate(variants):
        passed = detector.detect(variant)
        print("Variant", count, variant.getChromName(), variant.getPos(),
            "AF=", variant.getAF(), "passed=", passed)

    detector.close()

    with open(outfilename, "w") as outp:
        for variant in variants:
            if not variant.getProp("PASSED"):
                continue
            print("\t".join([variant.getChromName(), str(variant.getPos()),
                variant.getRef(), variant.getAlt(),
                "%.05f" % variant.getProp("PP"),
                "%.05f" % variant.getProp("AF_unrel")]), file = outp)


#========================================
if __name__ == "__main__":
    from argparse import ArgumentParser
    parser = ArgumentParser(
        description = "Detect DeNovo 2nd stage (high precision) algorithm")
    parser.add_argument("--output", "-O",
        help = "Output resulting file, required", required = True)
    parser.add_argument("--initial", "-I",
        help = "Initial file with 1st stage result, required",
        required = True)
    parser.add_argument("--unrelated", "-U",
        help = "Path to directory with libraries of unrelated, required")
    parser.add_argument("--trio", "-T",
        help = "Path to file describing trio BAM-files, required",
        required = True)
    parser.add_argument("--dump", "-D",
        help = "Tool for test original Anvoy code")
    parser.add_argument("--mdl", "-M",
        help = "AD-MDL index file, can be used instead of -U")

    args = parser.parse_args()

    if (args.unrelated,  args.mdl) == (None,  None):
        print("Either -I or -U option is required")

    runner(args.output, args.initial, args.unrelated, args.mdl,
        args.trio, args.dump)
