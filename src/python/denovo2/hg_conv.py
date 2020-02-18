#  Copyright (c) 2019. Partners HealthCare, Harvard Medical School’s
#  Department of Biomedical Informatics, Sergey Trifonov
#
#  Developed by Sergey Trifonov and Michael Bouzinier,
#  based on contributions by Anwoy Kumar Mohanty, Andrew Bjonnes,
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

import sys
from pyliftover import LiftOver

#========================================
class Hg19_38:
    sHandler = None
    sChromMap = dict()

    @classmethod
    def prepareChrom(cls, chrom):
        if chrom.startswith("chr"):
            return chrom
        chrom = str(chrom).upper()
        chrom = {"0": "M", "23": "X", "24": "Y"}.get(chrom, chrom)
        return "chr" + chrom

    @classmethod
    def convertPos(cls, chrom, pos):
        if cls.sHandler is None:
            print("Initializing hg19 -> hg38 liftover conversion",
                file = sys.stderr)
            cls.sHandler = LiftOver("hg19", "hg38")
        if chrom not in cls.sChromMap:
            cls.sChromMap[chrom] = cls.prepareChrom(chrom)
        try:
            coord  = cls.sHandler.convert_coordinate(cls.sChromMap[chrom], pos - 1)
        except Exception:
            return None
        if (len(coord) == 0):
            return None
        return coord[0][1] + 1

#========================================
class Hg38_19:
    sHandler = None
    sChromMap = dict()

    @classmethod
    def prepareChrom(cls, chrom):
        if chrom.startswith("chr"):
            return chrom
        chrom = str(chrom).upper()
        chrom = {"0": "M", "23": "X", "24": "Y"}.get(chrom, chrom)
        return "chr" + chrom

    @classmethod
    def convertPos(cls, chrom, pos):
        if cls.sHandler is None:
            print("Initializing hg38 -> hg19 liftover conversion",
                file = sys.stderr)
            cls.sHandler = LiftOver("hg38", "hg19")
        if chrom not in cls.sChromMap:
            cls.sChromMap[chrom] = cls.prepareChrom(chrom)
        try:
            coord  = cls.sHandler.convert_coordinate(cls.sChromMap[chrom], pos - 1)
        except Exception:
            return None
        if (len(coord) == 0):
            return None
        return coord[0][1] + 1
