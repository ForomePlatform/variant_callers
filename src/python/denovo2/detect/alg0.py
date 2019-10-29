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
# ALT_count = 1 !
# combos = (ALT_count + 1) * (ALT_count + 2) / 2 = 2 * 3 /2 =3
# GT_ordering_alternate(ALT_count) = np.array(
#    [[ 0., 0.], [ 0., 1.], [ 1.,  1.]])
# GT_likelihood_wrt_allele_calc(ALT_count = 1) =
#   np.array([[ 1.,  0. ], [ 0.5,  0.5], [ 0.,  1. ]])

GT_likelihood_L  = np.log(
    np.array([[1., 1E-100], [0.5, 0.5], [1E-100, 1.]]))
GT_likelihood_L_T = GT_likelihood_L.T

#========================================
# M1_L_calc
def prepRhoMatr(rho):
    side_v = np.log(1. - rho)
    main_v = np.log(rho)
    return [np.array([main_v, side_v]), np.array([side_v, main_v])]

#========================================
def logMean(v):
    return np.log(np.mean(np.exp(v)) + 1E-100)

def normLogVec(log_vec):
    vec = np.exp(np.array(log_vec).flatten() - np.max(log_vec))
    sum_val = vec.sum()
    return vec/sum_val, sum_val

#========================================
def makeMM2(MM1):
    global GT_likelihood_L
    # Anvoy's code: + np.log(2)
    return np.matrix([[logMean(MM1[i] + GT_likelihood_L[j])
        for j in range(3)] for i in range(2)])

def makeMM4(GT_marg, MM2):
    global GT_likelihood_L_T
    # Anvoy's code: + np.log(3)
    return np.array([[
        logMean(GT_likelihood_L_T[i] + GT_marg - MM2[j])
        for j in range(2)] for i in range(2)])

#========================================
def T_term_calc_for_rho(A_marg_L, AD):
    assert len(A_marg_L) == AD.size, "ERROR in T_term_calc"

    print("A_marg_L=", A_marg_L)
    T1_term, T2_term = 0., 0.
    for k, ad_k in enumerate(AD):
        A_marg_L_k = A_marg_L[k]
        A_marg_temp = np.exp(A_marg_L_k - np.max(A_marg_L_k))
        A_marg = A_marg_temp / np.sum(A_marg_temp)
        print("k=", k, "A_marg=", A_marg)
        T1_term += A_marg[k] * ad_k
        T2_term += (1. - A_marg[k]) * ad_k
    return T1_term, T2_term

#========================================
def EM_step(ADf_list, ADr_list, rho_f, rho_r, prior_L, allele_freq):

    assert len(ADf_list) == len(ADr_list), "ERROR1 in EM_step"

    AF = max(0., allele_freq)

    T1_f, T2_f = 1., 1.
    T1_r, T2_r = 1., 1.
    f0, f2 = (1. - AF) ** 2, AF ** 2
    T_for_prior = np.array([f0, 1. - f0 - f2, f2]) * 1000. + 1.
    joint_probty = (
        np.log(rho_f) + np.log(1. - rho_f) +
        np.log(rho_r) + np.log(1.- rho_r) +
        np.sum(T_for_prior * prior_L))

    MM1_f = prepRhoMatr(rho_f)
    MM2_f = makeMM2(MM1_f)

    MM1_r = prepRhoMatr(rho_r)
    MM2_r = makeMM2(MM1_r)

    for ADf, ADr in zip(ADf_list, ADr_list):
        GT_marg_L = prior_L + ADf * MM2_f + ADr * MM2_r

        MM4_f = makeMM4(GT_marg_L, MM2_f)
        MM4_r = makeMM4(GT_marg_L, MM2_r)

        marg_l_max = np.max(GT_marg_L)
        GT_marg = np.exp(np.array(GT_marg_L - marg_l_max).flatten())
        gt_sum = np.sum(GT_marg)

        T_for_prior += GT_marg / gt_sum
        joint_probty += np.log(gt_sum) + marg_l_max

        TT_f = T_term_calc_for_rho((MM1_f + MM4_f).T, ADf)
        T1_f += TT_f[0]
        T2_f += TT_f[1]
        TT_r = T_term_calc_for_rho((MM1_r + MM4_r).T, ADr)
        T1_r += TT_r[0]
        T2_r += TT_r[1]

    rho_f_new = 1. / (1. + T2_f / T1_f)
    rho_r_new = 1. / (1. + T2_r / T1_r)
    prior_new = T_for_prior / np.sum(T_for_prior)
    prior_L_new = np.log(prior_new)

    return rho_f_new, rho_r_new, prior_L_new, joint_probty

#========================================
def EM_full(ADfs, ADrs, rho_f, rho_r, prior_L, allele_freq):
    joint_probty_s = []
    while (len(joint_probty_s) < 3 or
            np.abs(joint_probty_s[-2] - joint_probty_s[-1]) > 1E-7):
        rho_f, rho_r, prior_L, joint_probty_new = \
            EM_step(ADfs, ADrs, rho_f, rho_r, prior_L, allele_freq)
        joint_probty_s.append(joint_probty_new)
    return rho_f, rho_r, prior_L

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

#========================================
def GTL_L_calc(ADf, ADr, rho_f, rho_r):
    MM1_f = prepRhoMatr(rho_f)
    MM2_f = makeMM2(MM1_f)

    MM1_r = prepRhoMatr(rho_r)
    MM2_r = makeMM2(MM1_r)

    GTL_L = np.array(ADf * MM2_f + ADr * MM2_r).flatten()
    return GTL_L - np.max(GTL_L)

#========================================
# denovo_P_calc
def evalResult(ADfs, ADrs, rho_f, rho_r, table_L, prior_L):
    M_GL_L, D_GL_L, C_GL_L = [GTL_L_calc(ADfs[idx] ,ADrs[idx],
        rho_f, rho_r) for idx in range(3)]
    combos = 3
    work_column = np.empty(combos ** 2)
    for idx1 in range(combos):
        for idx2 in range(combos):
            work_column[(idx1 * combos) + idx2] = (
                prior_L[idx1] + M_GL_L[idx1] +
                prior_L[idx2] +  + D_GL_L[idx2])

    work_table = (table_L + np.tile(C_GL_L,[combos ** 2, 1]) +
        np.tile(np.reshape(work_column, [combos**2, 1]),[1, combos]))
    work_table -= np.max(work_table)
    work_table = np.exp(work_table)
    work_table /= np.sum(work_table)
    PP = np.max(np.array([work_table[0][1], work_table[0][2]]))
    return PP

#========================================
#PP_calc
def calculateVariant(variant,
    trio_samfiles, unrelated_samfiles, debug_mode):

    ADfs_U, ADrs_U = unrelated_samfiles.mineAD(variant)

    if ADfs_U is None:
        variant.setProp("PASSED", False)
        return

    if debug_mode:
        variant.setProp("ADfs_U", ADfs_U)
        variant.setProp("ADrs_U", ADrs_U)

    rho_f_old = 0.8
    rho_r_old = 0.8
    prior_L_old = np.log(np.array([1./3, 1./3, 1./3]))

    rho_f_new, rho_r_new, prior_L_new = EM_full(
        ADfs_U, ADrs_U, rho_f_old, rho_r_old, prior_L_old,
        variant.getProp("AF"))
    variant.setProp("rho_f", rho_f_new)
    variant.setProp("rho_r", rho_r_new)
    variant.setProp("prior_L", prior_L_new)

    print("Model:", rho_f_new, rho_r_new, prior_L_new)

    AF_unrel = 0.
    for i in range(ADfs_U.shape[0]):
        temp1 = GTL_L_calc(ADfs_U[i], ADrs_U[i], rho_f_new, rho_r_new)
        temp = temp1 + prior_L_new
        temp -= np.max(temp)
        temp = np.exp(temp)
        temp /= np.sum(temp)
        AF_unrel += temp[1] + temp[2] * 2.
    AF_unrel /= 2. * ADfs_U.shape[0]
    variant.setProp("AF_unrel", AF_unrel)

    print("AF_unrel=", AF_unrel, "AF=", variant.getProp("AF"))

    ADfs, ADrs = trio_samfiles.mineAD(variant)
    assert len(ADfs) == 3 and len(ADrs) == 3
    if debug_mode:
        variant.setProp("ADfs", ADfs)
        variant.setProp("ADrs", ADrs)

    #PASS:if AF_unrel<0.01 and ALT_read_checker_in_parents(ADfs,ADrs):
    if (AF_unrel < 0.01 and
            ADfs[0][1] + ADfs[1][1] + ADrs[0][1] + ADrs[1][1] <= 3):
        variant.setProp("PASSED", True)
    else:
        variant.setProp("PASSED", False)
        print ("Not passed", ADfs[0][1] + ADfs[1][1] + ADrs[0][1] + ADrs[1][1]);
        return

    table_L = np.log(makeTableGen(1e-8))
    PP = evalResult(ADfs, ADrs, rho_f_new, rho_r_new, table_L, prior_L_new)
    variant.setProp("PP", PP)
