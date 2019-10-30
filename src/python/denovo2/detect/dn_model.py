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

import numpy as np
from .algorithm import RhoModel, EM_full, evalAF, evalPP, normLogMatr
from adlib.mdl_io import MDL_Reader

#========================================
class DeNovo_Model:
    sStartPriorL = np.log(np.array([1./3, 1./3, 1./3]))
    sVerboseMode = False

    def __init__(self, rho_model, prior_L, af_unrel = None):
        self.mRhoModel = rho_model
        self.mPriorL = prior_L
        self.mAF_unrel = af_unrel

    def isBad(self):
        return self.mRhoModel is None

    def evalVariant(self, variant, trio_samfiles, report_mode = True):
        if self.mRhoModel is None:
            variant.setProp("PASSED", False)
            return

        ADfs, ADrs = trio_samfiles.mineAD(variant)
        assert len(ADfs) == 3 and len(ADrs) == 3

        #PASS:if AF_unrel<0.01 and ALT_read_checker_in_parents(ADfs,ADrs):
        if (ADfs[0][1] + ADfs[1][1] + ADrs[0][1] + ADrs[1][1] > 3):
            variant.setProp("PASSED", False)
            if self.sVerboseMode:
                print ("Not passed by trio data:", ADfs, ADrs,
                    ADfs[0][1] + ADrs[0][1], ADfs[1][1] + ADrs[1][1])
            return

        variant.setProp("PASSED", True)
        if report_mode is True:
            rho_f, rho_r = self.mRhoModel.getRhoPair()
            variant.setProp("rho_f", rho_f)
            variant.setProp("rho_r", rho_r)
            variant.setProp("prior_L", self.mPriorL)
            if self.mAF_unrel is not None:
                variant.setProp("AF_unrel", self.mAF_unrel)
        PP = evalPP(ADfs, ADrs, self.mRhoModel, self.mPriorL)
        variant.setProp("PP", PP)

    @classmethod
    def createByADLib(cls, variant, ad_lib):
        ADfs_U, ADrs_U = ad_lib.mineAD(variant)

        if ADfs_U is None:
            return DeNovo_Model(None, None, None)

        rho_mod = RhoModel.create(ADfs_U, ADrs_U)
        prior_L = EM_full(ADfs_U, ADrs_U, rho_mod,
            cls.sStartPriorL, variant.getAF())

        AF_unrel = evalAF(ADfs_U, ADrs_U, rho_mod, prior_L)
        if AF_unrel >= 0.01:
            if cls.sVerboseMode:
                print("Not passed by AF_unrel: %.05f" % AF_unrel)
            return DeNovo_Model(None, None, None)
        return DeNovo_Model(rho_mod, prior_L, AF_unrel)

    @classmethod
    def makeApproxPosModel(cls, ADfs, ADrs):
        if (ADfs is None or len(ADfs) < 1 or
                (max([ad.max() for ad in ADfs + ADrs])) == 0):
            return [0]* 8
        rho_mod = RhoModel.create(ADfs, ADrs)
        pos_model = rho_mod.getRhoPair()
        for af in (0., .05):
            prior_L = EM_full(ADfs, ADrs, rho_mod,
                DeNovo_Model.sStartPriorL, af)
            prior, _ = normLogMatr(prior_L)
            pos_model += [prior[0], prior[1],
                evalAF(ADfs, ADrs, rho_mod, prior_L)]
        return [prob2ushort(val) for val in pos_model]

    @classmethod
    def createByApproxPosModel(cls, variant, pos_ushort_model):
        if pos_ushort_model is None:
            if cls.sVerboseMode:
                print("Not passed - no data")
            return DeNovo_Model(None, None, None)

        rho_f, rho_r, p0_0, p2_0, afu_0, p0_1, p2_1, afu_1 = [
            ushort2prob(val) for val in pos_ushort_model]
        af_cur = max(0., variant.getAF())
        if af_cur > .05:
            if cls.sVerboseMode:
                print("Not passed by AF: %.05f" % af_cur)
            return DeNovo_Model(None, None, None)

        q1 = af_cur / .05
        q0 = 1. - q1
        afu = q0 * afu_0 + q1 * afu_1
        if afu >= .01:
            if cls.sVerboseMode:
                print("Not passed by AF_unrel: %.05f" % afu)
            return DeNovo_Model(None, None, None)

        p0  = q0 * p0_0  + q1 * p0_1
        p2  = q0 * p2_0  + q1 * p2_1
        prior_L = np.log(np.array([p0, 1 - p0 - p2, p2]))
        return DeNovo_Model(RhoModel(rho_f, rho_r), prior_L, afu)

#========================================
class DeNovo_MDL_Reader(MDL_Reader):
    PREFIX = "AD/MDL.v1.0"

    def __init__(self, fname):
        MDL_Reader.__init__(self, fname, self.PREFIX, 'H', 8)

    def getPosModel(self, variant):
        return DeNovo_Model.createByApproxPosModel(variant,
            self.getModel(variant.getChromNum(), variant.getPos()))

#========================================
MAX_USHORT = 50000
def prob2ushort(val):
    global MAX_USHORT
    return min(MAX_USHORT, max(0, int(val * MAX_USHORT)))

def ushort2prob(val):
    global MAX_USHORT
    return min(1., max(0., float(val) / MAX_USHORT))

#========================================
