#  Copyright (c) 2019. Sergey Trifonov, Michael Bouzinier, Anwoy Kumar Mohanty,
#  Ignat Leshiner, Shamil Sunyaev
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
#========================================
def normLogMatr(log_matr):
    log_max = np.max(log_matr)
    vec = np.exp(np.array(log_matr).flatten() - log_max)
    sum_val = vec.sum()
    return vec/sum_val, log_max + np.log(sum_val)

#========================================
#========================================
def evalRho(AD_list):
    TT = np.ones(2)
    for AD in AD_list:
        TT += AD
    return TT[0]/TT.sum()

#========================================
def rho2prior(rho):
    return np.log(np.matrix([
        [rho,  0.5 ,  1. - rho], [1. - rho,  0.5 ,  rho]]))

#========================================
class RhoModel:
    def __init__(self, rho_f, rho_r):
        self.mRhoF = rho_f
        self.mRhoR = rho_r
        self.mPriorF = rho2prior(self.mRhoF)
        self.mPriorR = rho2prior(self.mRhoR)

    def apply(self, ADf, ADr):
        return np.array(ADf * self.mPriorF + ADr * self.mPriorR).flatten()

    def getRhoPair(self):
        return [self.mRhoF, self.mRhoR]

    @staticmethod
    def create(ADf_list, ADr_list):
        return RhoModel(evalRho(ADf_list), evalRho(ADr_list))

#========================================
#========================================
def EM_step(ADf_list, ADr_list, rho_mod, prior_L, allele_freq):
    AF = max(0., allele_freq)
    f0, f2 = (1. - AF) ** 2, AF ** 2
    T_for_prior = np.array([f0, 1. - f0 - f2, f2]) * 1000. + 1.
    joint_probty = np.sum(T_for_prior * prior_L)

    for ADf, ADr in zip(ADf_list, ADr_list):
        gt_marg, gt_log = normLogMatr(prior_L + rho_mod.apply(ADf, ADr))
        T_for_prior += gt_marg
        joint_probty += gt_log

    prior_L_new = np.log(T_for_prior / np.sum(T_for_prior))
    return prior_L_new, joint_probty

#========================================
def EM_full(ADfs, ADrs, rho_mod, prior_L, allele_freq):
    joint_probty_seq = []
    while (len(joint_probty_seq) < 3 or
            np.abs(joint_probty_seq[-2] - joint_probty_seq[-1]) > 1E-7):
        prior_L, joint_probty = EM_step(ADfs, ADrs, rho_mod,
            prior_L, allele_freq)
        joint_probty_seq.append(joint_probty)
    return prior_L

#========================================
# table_gen(1, 1e-8)
def makeTableGen3(dd = 1e-8):
    return np.array([
        [[  1. - (2 * dd),   (2 - dd) * dd,   dd * dd],
        [  .5 - (dd / 2),   .5,   dd/2],
        [  dd * (1. - dd),   1. - 2 * dd * (1. - dd), dd * (1.- dd)]],
        [[  .5 - (dd / 2),   .5,   dd/2],
        [  .25,   .5,   .25],
        [  dd/2,   .5,   .5 - (dd / 2)]],
        [[  dd * (1 - dd),   1. - 2 * dd * (1. - dd), dd * (1 - dd)],
        [  dd/2,   .5,   .5 - (dd / 2)],
        [  dd * dd,   (2 - dd) * dd,   1. - (2 * dd)]]])

sTableGen_L3 = np.log(makeTableGen3())

#========================================
# denovo_P_calc
def evalPP(ADfs, ADrs, rho_mod, prior_L):
    fam_gtl = [rho_mod.apply(ADf, ADr)
        for ADf, ADr in zip(ADfs, ADrs)]

    work_table = sTableGen_L3.copy()
    for i in range(3):
        for j in range(3):
            for k in range(3):
                work_table[i, j, k] += (fam_gtl[2][k] +
                    prior_L[j] + fam_gtl[1][j] +
                    prior_L[i] + fam_gtl[0][i])
    norm_vec, _ = normLogMatr(work_table)
    return max(norm_vec[1], norm_vec[2])

#========================================
def evalAF(ADfs, ADrs, rho_mod, prior_L):
    af_sum = 0.
    for ADf, ADr in zip(ADfs, ADrs):
        prob_vec, _ = normLogMatr(prior_L + rho_mod.apply(ADf, ADr))
        af_sum += prob_vec[1] + 2 * prob_vec[2]
    return af_sum / (2. * ADfs.shape[0])

