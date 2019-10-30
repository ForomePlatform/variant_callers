#  Copyright (c) 2019. Sergey Trifonov, Michael Bouzinier
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
from datetime import datetime
from .read_pysam import AD_LibCollection
from adlib.mdl_io import MDL_Builder
from .dn_model import DeNovo_Model, DeNovo_MDL_Reader
from .detect2 import VariantHandler

#========================================
def buildApproxModel(fname, ad_lib_coll, report_count = -1):
    mdl_builder = MDL_Builder(fname, DeNovo_MDL_Reader.PREFIX, 'H', 8)
    portion_count = 0
    dt0 = datetime.now()
    while True:
        portion_info = ad_lib_coll._nextPortions()
        if portion_info is None:
            break
        portion_count += 1
        if report_count > 0 and portion_count % report_count == 0:
            dt1 = datetime.now()
            print("portions=", portion_count, " time=", dt1 - dt0,
                file = sys.stderr)
            dt0 = dt1
        chrom, pos0, pos1 = portion_info
        for pos in range(pos0, pos1):
            ADfs_U, ADrs_U = ad_lib_coll.mineAD(VariantHandler(
                None, pos, None, None, .0, chrom_num = chrom))
            pos_model = DeNovo_Model.makeApproxPosModel(ADfs_U, ADrs_U)
            mdl_builder.addRecord(chrom, pos, pos_model)
    mdl_builder.close()
    if report_count > 0:
        print("Total %d portions" % portion_count, file = sys.stderr)

def debugBuild(ad_lib_dir, out_mdl_name):
    variants = [
        VariantHandler("chr1",  12837210, "C", "T", 0.000329),
        VariantHandler("chr15", 83013622, "A", "G", 0.000283),
        VariantHandler("chr16", 75536430, "C", "T", -1.0),
        VariantHandler("chr17", 39633973, "G", "A", -1.0)]

    class _DebugLibReader(AD_LibCollection):
        def __init__(self, dirname, variants):
            AD_LibCollection.__init__(self, dirname)
            self.mVariants = variants[:]

        def _nextPortions(self):
            if len(self.mVariants) == 0:
                return None
            var = self.mVariants.pop(0)
            chrom, pos = var.getChromNum(), var.getPos()
            return (chrom, pos - 100, pos + 100)

    ad_lib = _DebugLibReader(ad_lib_dir, variants)
    buildApproxModel(out_mdl_name, ad_lib, 1)

#========================================
if __name__=="__main__":
    if sys.argv[1] == "--test":
        ad_lib_dir, out_mdl_name = sys.argv[2:]
        debugBuild(ad_lib_dir, out_mdl_name)
    else:
        ad_lib_dir, out_mdl_name = sys.argv[1:]
        ad_lib = AD_LibCollection(ad_lib_dir)
        buildApproxModel(out_mdl_name, ad_lib, 1)



