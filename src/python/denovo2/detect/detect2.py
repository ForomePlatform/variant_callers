import numpy as np

from .algorithm import calculateVariant
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
        chrom_num = {"M":0, "MT": 0, "X":23, "Y":24}.get(chrom_name, 25)
    assert chrom_num < 25, "Bad chromosome: " + chrom_name
    return chrom_num

def num2chrom(chrom_num):
    return (str(chrom_num) if 0 < chrom_num <= 22
        else {0: "MT", 23: "X", 24: "Y"}[chrom_num])

#========================================
# Variant works
#========================================
class VariantHandler:
    def __init__(self, line):
        self.mLine = line.strip()
        params = self.mLine.split()
        self.mChromName = params[0]
        self.mChromNum = chrom2num(self.mChromName)
        self.mPos = int(params[1])
        self.mRef = params[3]
        self.mAlt = params[4]
        self.mProps = dict()
        for pair in params[7].split(';'):
            arg, _, val = pair.partition('=')
            if arg == "ExAC_AF_computed":
                self.mProps["AF"] = float(val)

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

    def getLine(self):
        return self.mLine

    def getProp(self, key):
        return self.mProps.get(key)

    def setProp(self, key, val):
        self.mProps[key] = val

    @classmethod
    def loadFile(cls, filename):
        variants = []
        with open(filename, "r") as inp:
            for line in inp:
                var = VariantHandler(line)
                if var.isGood():
                    variants.append(var)
        variants.sort(key = lambda var: (var.getChromNum(), var.getPos()))
        return variants

#========================================
def runner(outfilename, initial_filename, unrelated_dir,
        trio_filename, dump_filename):
    variants = VariantHandler.loadFile(initial_filename)
    unrelated_samfiles = AD_LibCollection(unrelated_dir, dump_filename)
    trio_samfiles = PysamList(trio_filename)

    for count, variant in enumerate(variants):
        print("Variant", count, variant.getChromName(), variant.getPos(),
            "AF=", variant.getProp("AF"))
        if dump_filename:
            unrelated_samfiles.mineAD(variant)
        calculateVariant(variant,
            trio_samfiles, unrelated_samfiles)
        if variant.getProp("PASSED"):
            print("passed")
    unrelated_samfiles.finishUp()

    cnt = 0
    with open(outfilename, "w") as outp:
        for variant in variants:
            if not variant.getProp("PASSED"):
                continue
            cnt += 1
            print("\t".join([str(cnt), variant.getChromName(),
                str(variant.getPos()), variant.getRef(), variant.getAlt(),
                "AF=%r" % variant.getProp("AF"),
                "rhos=%r,%r" % (variant.getProp("rho_f"),
                    variant.getProp("rho_r")),
                "prior=%r" % np.exp(variant.getProp("prior_L")),
                "PP=%r" % variant.getProp("PP"),
                "AF_unrel=%r" % variant.getProp("AF_unrel"),
                "CSQ_gene=%r" % variant.getProp("CSQ_gene")]), file = outp)

#========================================
if __name__=="__main__":
    from argparse import ArgumentParser
    parser = ArgumentParser(
        description = "Detect DeNovo 2nd stage (high precision) algorithm")
    parser.add_argument("--output", "-O",
        help = "Output resulting file, required", required=True)
    parser.add_argument("--initial", "-I",
        help = "Initial file with 1st stage result, required", required=True)
    parser.add_argument("--unrelated", "-U",
        help = "Path to directory with libraries of unrelated, required", required=True)
    parser.add_argument("--trio", "-T",
        help = "Path to file describing trio BAM-files, required", required=True)
    parser.add_argument("--dump", "-D",
        help = "Tool for test original Anvoy code")

    args = parser.parse_args()

    print("outfilename=", args.output)
    print("initial_filename=", args.initial)
    print("unrelated_dir=", args.unrelated)
    print("trio_filename=", args.trio)

    runner(args.output, args.initial, args.unrelated, args.trio,
        args.dump)
