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
GT_likelihood_L  = np.log(
    np.array([[1., 1E-100], [0.5, 0.5], [1E-100, 1.]]))

#========================================
def logMean(v):
    return np.log(np.mean(np.exp(v)) + 1E-100)

def normLogMatr(log_matr):
    log_max = np.max(log_matr)
    vec = np.exp(np.array(log_matr).flatten() - log_max)
    sum_val = vec.sum()
    return vec/sum_val, log_max + np.log(sum_val)

#========================================
#========================================
def makeMM2(rho):
    global GT_likelihood_L
    main_v = np.log(rho)
    side_v = np.log(1. - rho)
    MM1 = [np.array([main_v, side_v]), np.array([side_v, main_v])]

    return np.matrix([[logMean(MM1[i] + GT_likelihood_L[j])
        for j in range(3)] for i in range(2)])

#========================================
def eval_rho(AD_list):
    TT = np.ones(2)
    for AD in AD_list:
        TT += AD
    return 1. / (1. + TT[1] / TT[0])

#========================================
def EM_step(ADf_list, ADr_list, rho_f, rho_r, prior_L, allele_freq):

    assert len(ADf_list) == len(ADr_list), "ERROR1 in EM_step"

    AF = max(0., allele_freq)

    f0, f2 = (1. - AF) ** 2, AF ** 2
    T_for_prior = np.array([f0, 1. - f0 - f2, f2]) * 1000. + 1.
    joint_probty = (
        np.log(rho_f) + np.log(1. - rho_f) +
        np.log(rho_r) + np.log(1.- rho_r) +
        np.sum(T_for_prior * prior_L))

    MM2_f = makeMM2(rho_f)
    MM2_r = makeMM2(rho_r)

    for ADf, ADr in zip(ADf_list, ADr_list):
        GT_marg_L = prior_L + ADf * MM2_f + ADr * MM2_r

        GT_marg, gt_log = normLogMatr(GT_marg_L)
        T_for_prior += GT_marg
        joint_probty += gt_log

    prior_L_new = np.log(T_for_prior / np.sum(T_for_prior))

    return prior_L_new, joint_probty

#========================================
def EM_full(ADfs, ADrs, rho_f, rho_r, prior_L, allele_freq):
    joint_probty_s = []
    while (len(joint_probty_s) < 3 or
            np.abs(joint_probty_s[-2] - joint_probty_s[-1]) > 1E-7):
        prior_L, joint_probty_new = EM_step(ADfs, ADrs, 
            rho_f, rho_r, prior_L, allele_freq)
        joint_probty_s.append(joint_probty_new)
    return prior_L

#========================================
def GTL_L_calc(ADf, ADr, rho_f, rho_r):
    MM2_f = makeMM2(rho_f)
    MM2_r = makeMM2(rho_r)
    GTL_L = np.array(ADf * MM2_f + ADr * MM2_r).flatten()
    return GTL_L - np.max(GTL_L)

#========================================
# table_gen(1, 1e-8)
def makeTableGen(dd = 1e-8):
    return np.array([
        [  1. - (2 * dd),   (2 - dd) * dd,   dd * dd],
        [  .5 - (dd / 2),   .5,   dd/2],
        [  dd * (1. - dd),   1. - 2 * dd * (1. - dd), dd * (1.- dd)],
        [  .5 - (dd / 2),   .5,   dd/2],
        [  .25,   .5,   .25],
        [  dd/2,   .5,   .5 - (dd / 2)],
        [  dd * (1 - dd),   1. - 2 * dd * (1. - dd), dd * (1 - dd)],
        [  dd/2,   .5,   .5 - (dd / 2)],
        [  dd * dd,   (2 - dd) * dd,   1. - (2 * dd)]])

sTableGen_L1 = np.log(makeTableGen())

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
def evalResult(ADfs, ADrs, rho_f, rho_r, prior_L):
    fam_gtl = [GTL_L_calc(ADfs[idx], ADrs[idx], rho_f, rho_r) 
        for idx in range(3)]
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
#PP_calc
def calculateVariant(variant, trio_samfiles, unrelated_samfiles):

    ADfs_U, ADrs_U = unrelated_samfiles.mineAD(variant)

    if ADfs_U is None:
        variant.setProp("PASSED", False)
        return

    rho_f = eval_rho(ADfs_U)
    rho_r = eval_rho(ADrs_U)
    variant.setProp("rho_f", rho_f)
    variant.setProp("rho_r", rho_r)
    
    prior_L_start = np.log(np.array([1./3, 1./3, 1./3]))

    prior_L = EM_full(ADfs_U, ADrs_U, rho_f, rho_r, prior_L_start,
        variant.getProp("AF"))
    variant.setProp("prior_L", prior_L)

    AF_unrel = 0.
    for i in range(ADfs_U.shape[0]):
        prob_vec, _ = normLogMatr(prior_L +
            GTL_L_calc(ADfs_U[i], ADrs_U[i], rho_f, rho_r))
        AF_unrel += prob_vec[1] + 2 * prob_vec[2]
    AF_unrel /= 2. * ADfs_U.shape[0]
    variant.setProp("AF_unrel", AF_unrel)

    ADfs, ADrs = trio_samfiles.mineAD(variant)
    assert len(ADfs) == 3 and len(ADrs) == 3

    #PASS:if AF_unrel<0.01 and ALT_read_checker_in_parents(ADfs,ADrs):
    if (AF_unrel < 0.01 and
            ADfs[0][1] + ADfs[1][1] + ADrs[0][1] + ADrs[1][1] <= 3):
        variant.setProp("PASSED", True)
    else:
        variant.setProp("PASSED", False)
        print ("Not passed", ADfs[0][1] + ADfs[1][1] + ADrs[0][1] + ADrs[1][1]);
        return

    PP = evalResult(ADfs, ADrs, rho_f, rho_r, prior_L)
    variant.setProp("PP", PP)
