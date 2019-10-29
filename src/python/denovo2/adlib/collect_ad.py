import sys, pysam, time
from .hg19_rd import Hg19_Reader
from .hg19_portion import Hg19_Portion
from .ad_person import AD_PersonData
from .ad_miner import AD_PortionMiner, AD_PersonDataWriter

#========================================
if len(sys.argv) < 4:
    print("More arguments required", file = sys.stderr)
    sys.exit(1)

if sys.argv[1] == "-donothing":
    print(sys.argv, file = sys.stderr)
    sys.exit(0)

rd = Hg19_Reader(sys.argv[1], chrom_is_int = True, upper_case = True)
samfile = pysam.AlignmentFile(sys.argv[2], "rb")
writer = AD_PersonDataWriter(sys.argv[3])

if len(sys.argv) > 4 and not "--" in sys.argv[4]:
    print("Merge with:", sys.argv[4], file = sys.stderr)
    merge_reader = AD_PersonData(sys.argv[4], True)
else:
    merge_reader = None

cnt = 0
while rd.read():
    cnt += 1
    if cnt % 10 == 0:
        print("Progress:", cnt, "portions...", file = sys.stderr)
    tm0 = time.time()
    hg19_portion = Hg19_Portion(rd.getCurChrom(),
        rd.getCurDiap(), rd.getCurLetters())
    ad_portion = None
    if merge_reader is not None:
        ad_portion = merge_reader.directReadPortion()
        if ad_portion is False:
            print("Merge ad-data portion failed", file = sys.stderr)
            ad_portion = None
        elif ad_portion is None:
            print("Merge ad-data ends", file = sys.stderr)
            merge_reader = None
    if ad_portion is None:
        ad_portion = AD_PortionMiner(hg19_portion, samfile)
    else:
        tm0 = None
    writer.addPortion(ad_portion)
    if tm0 is not None:
        print(ad_portion.report(), ("sec: %0.1f" % (time.time() - tm0)),
            file = sys.stderr)
writer.close()

print("Done:", cnt, file = sys.stderr)

