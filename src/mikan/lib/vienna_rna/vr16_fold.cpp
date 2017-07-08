#include <iostream>
#include <algorithm>
#include "vr16_energy.hpp"
#include "vr16_fold.hpp"

namespace vr16 {

//
// VR16Fold methods
//
int VR16Fold::init_arrays(int pStrLen) {
    if (pStrLen < 1) {
        std::cerr << "Error: Length must be greater 0" << std::endl;
        return 1;
    }

    if (mInitLength > 0) {
        clear_arrays();
    }

    resize_arrays(pStrLen);

    mIdex.resize((unsigned) pStrLen + 1);
    for (int i = 1; i <= pStrLen; ++i) {
        mIdex[i] = (i * (i - 1)) >> 1; /* n(n-1)/2 */
    }

    mInitLength = pStrLen;

    return 0;
}

void VR16Fold::clear_arrays() {
    mF5.clear();
    mPtype.clear();

    mOpts.mBasePair.clear();

    mInitLength = 0;
}

void VR16Fold::resize_arrays(int pStrLen) {
    mF5.resize((unsigned) pStrLen + 2);
    mPtype.resize((unsigned) (pStrLen * (pStrLen + 1)) / 2 + 2);

    mOpts.mBasePair.resize(1 + (unsigned) pStrLen / 2);
}

float VR16Fold::fold(std::string &pString, std::string &pStructure) {
    int energy;
    float energy2;
    int strLen;

    strLen = (int) pString.size();

    if (strLen > mInitLength) {
        init_arrays(strLen);
    }

    mBP.resize((unsigned) strLen + 2, 0);
    encode_seq(pString);
    make_ptypes(pStructure, mPtype, mBP, mOpts);

    energy = fill_arrays(pString);

    if (mOpts.mDoBacktrack) {
        mBT.backtrack(pString, 0, mPtype, mS1, mBP, mArrayC, mArrayML, mF5);
        mBT.parenthesis_structure(pStructure, strLen);
    }

    energy2 = (float) finalize_energy_calc(strLen, pStructure, energy, mOpts, mBP);

    mBP.clear();

    return energy2;

}

void VR16Fold::encode_seq(std::string &pString) {
    unsigned l;

    l = (unsigned) pString.size();

    mS.resize(l + 2);
    mS1.resize(l + 2);
    /* S1 exists only for the special X K and I bases and energy_set!=0 */
    mS[0] = l;
    mS1[0] = l;

    for (unsigned i = 1; i <= l; ++i) { /* make numerical encoding of sequence */
        mS[i] = mPairMat.encode_char(pString[i - 1], mOpts.mEnergySet);
        mS1[i] = mPairMat.get_alias(mS[i]); /* for mismatches of nostandard bases */
    }
    /* for circular folding add first base at position n+1 */
    mS[l + 1] = mS[1];
    mS1[l + 1] = mS1[1];

}

void VR16Fold::make_ptypes(
        std::string &pStructure,
        std::vector<char> &pPtype,
        std::vector<int> &pBP,
        VR16FoldOptions &pOpts) {
    int type_1 = 0;
    int ntype;
    int otype;
    int i;
    int j;

    int n = mS[0];
    for (int k = 1; k < n - TURN; ++k) {
        for (int l = 1; l <= 2; ++l) {
            ntype = 0;
            otype = 0;
            i = k;
            j = i + TURN + l;
            if (j > n) {
                continue;
            }
            type_1 = mPairMat.mPair[mS[i]][mS[j]];
            while ((i >= 1) && (j <= n)) {
                if ((i > 1) && (j < n)) {
                    ntype = mPairMat.mPair[mS[i - 1]][mS[j + 1]];
                }

                if (pOpts.mNoLonelyPairs && (!otype) && (!ntype)) {
                    type_1 = 0; /* i.j can only form isolated pairs */
                }

                pPtype[mIdex[j] + i] = (char) type_1;
                otype = type_1;
                type_1 = ntype;
                i--;
                j++;
            }
        }
    }

    if (pOpts.mFoldConstrained && (pStructure.size() != 0)) {
        int hx;
        std::vector<int> stack;
        char type;
        stack.resize((unsigned) n + 1);

        for (hx = 0, j = 1; j <= n; j++) {
            switch (pStructure[j - 1]) {
                case '|':
                    pBP[j] = -1;
                    break;
                case 'x': /* can't pair */
                    for (int l = 1; l < j - TURN; l++) {
                        pPtype[mIdex[j] + l] = 0;
                    }
                    for (int l = j + TURN + 1; l <= n; l++) {
                        pPtype[mIdex[l] + j] = 0;
                    }
                    break;
                case '(':
                    stack[hx++] = j;
                    for (int l = 1; l < j - TURN; l++) {
                        pPtype[mIdex[j] + l] = 0;
                    }
                    break;
                case '<': /* pairs upstream */
                    for (int l = 1; l < j - TURN; l++) {
                        pPtype[mIdex[j] + l] = 0;
                    }
                    break;
                case ')':
                    if (hx <= 0) {
                        std::cerr << "Error: Unbalanced brackets in constraints." << std::endl;
                        std::cerr << pStructure << std::endl;
                    }
                    i = stack[--hx];
                    type = pPtype[mIdex[j] + i];
                    for (int k = i + 1; k <= n; ++k) {
                        pPtype[mIdex[k] + i] = 0;
                    }
                    /* don't allow pairs i<k<j<l */
                    for (int l = j; l <= n; l++) {
                        for (int k = i + 1; k <= j; ++k) {
                            pPtype[mIdex[l] + k] = 0;
                        }
                    }
                    /* don't allow pairs k<i<l<j */
                    for (int l = i; l <= j; l++) {
                        for (int k = 1; k <= i; k++) {
                            pPtype[mIdex[l] + k] = 0;
                        }
                    }
                    for (int k = 1; k < j; k++) {
                        pPtype[mIdex[j] + k] = 0;
                    }
                    pPtype[mIdex[j] + i] = (type == 0) ? (char) 7 : (char) type_1;
                    for (int l = j + TURN + 1; l <= n; ++l) {
                        pPtype[mIdex[l] + j] = 0;
                    }
                    break;
                case '>': /* pairs downstream */
                    for (int l = j + TURN + 1; l <= n; ++l) {
                        pPtype[mIdex[l] + j] = 0;
                    }
                    break;
                default: // do nothing;
                    break;
            }
        }
        if (hx != 0) {
            std::cerr << "Error: Unbalanced brackets in constraints." << std::endl;
            std::cerr << pStructure << std::endl;
        }
        stack.clear();
    }
}

int VR16Fold::fill_arrays(std::string &pString) {
    /* fill "c", "fML" and "f5" arrays and return  optimal energy */
    int strLen = (int) pString.size();
    mArrayC.init_arrays(strLen);
    mArrayML.init_arrays(strLen);

    for (int i = strLen - TURN - 1; i >= 1; --i) { /* i,j in [1..length] */

        for (int j = i + TURN + 1; j <= strLen; ++j) {
            mArrayC.fill_arrays(i, j, mBP, mPtype, mS1, mArrayML.mFML, pString);
            /* done with c[i,j], now compute fML[i,j] */
            /* free ends ? -----------------------------------------*/
            mArrayML.fill_arrays(i, j, strLen, mPtype, mS1, mArrayC);
        }

        mArrayC.rotate_aux_idexes(strLen);
        mArrayML.reset_fmi(strLen);
    }

    calc_energy_of_fragments(strLen, mPtype, mF5, mOpts);

    return mF5[strLen];
}

void VR16Fold::calc_energy_of_fragments(
        int pStrLen,
        std::vector<char> &pPtype,
        std::vector<int> &pF5,
        VR16FoldOptions &pOpts) {
    int energy;
    int energy_t;
    int type_1;

    /* calculate energies of 5' and 3' fragments */
    pF5[TURN + 1] = 0;
    for (int j = TURN + 2; j <= pStrLen; ++j) {
        pF5[j] = pF5[j - 1];
        type_1 = pPtype[mIdex[j] + 1];
        if (type_1) {
            energy = mArrayC.mC[mIdex[j] + 1];
            if (type_1 > 2) {
                energy += mParams.mTerminalAU;
            }
            if ((pOpts.mDangles == 2) && (j < pStrLen)) /* double dangles */
            {
                energy += mParams.mDangle3[type_1][mS1[j + 1]];
            }
            if (energy < pF5[j]) {
                pF5[j] = energy;
            }
        }
        type_1 = pPtype[mIdex[j - 1] + 1];
        if ((type_1) && (pOpts.mDangles % 2 == 1)) {
            energy = mArrayC.mC[mIdex[j - 1] + 1] + mParams.mDangle3[type_1][mS1[j]];
            if (type_1 > 2) {
                energy += mParams.mTerminalAU;
            }
            if (energy < pF5[j]) {
                pF5[j] = energy;
            }
        }
        for (int i = j - TURN - 1; i > 1; --i) {
            type_1 = pPtype[mIdex[j] + i];
            if (type_1) {
                energy = pF5[i - 1] + mArrayC.mC[mIdex[j] + i];
                if (type_1 > 2) {
                    energy += mParams.mTerminalAU;
                }
                if (pOpts.mDangles == 2) {
                    energy += mParams.mDangle5[type_1][mS1[i - 1]];
                    if (j < pStrLen) {
                        energy += mParams.mDangle3[type_1][mS1[j + 1]];
                    }
                }
                if (energy < pF5[j]) {
                    pF5[j] = energy;
                }
                if (pOpts.mDangles % 2 == 1) {
                    energy = pF5[i - 2] + mArrayC.mC[mIdex[j] + i] + mParams.mDangle5[type_1][mS1[i - 1]];
                    if (type_1 > 2) {
                        energy += mParams.mTerminalAU;
                    }
                    if (energy < pF5[j]) {
                        pF5[j] = energy;
                    }
                }
            }
            type_1 = pPtype[mIdex[j - 1] + i];
            if ((type_1) && (pOpts.mDangles % 2 == 1)) {
                energy = mArrayC.mC[mIdex[j - 1] + i] + mParams.mDangle3[type_1][mS1[j]];
                if (type_1 > 2) {
                    energy += mParams.mTerminalAU;
                }
                energy_t = pF5[i - 1] + energy;
                if (energy_t < pF5[j]) {
                    pF5[j] = energy_t;
                }
                energy_t = pF5[i - 2] + energy + mParams.mDangle5[type_1][mS1[i - 1]];
                if (energy_t < pF5[j]) {
                    pF5[j] = energy_t;
                }
            }
        }
    }
}

double VR16Fold::finalize_energy_calc(
        int pStrLen,
        std::string &pStructure,
        int pEnergy,
        VR16FoldOptions &pOpts,
        std::vector<int> &pBP) {
    int bonus = 0;
    int bonus_cnt = 0;

    for (int i = 1; i <= pStrLen; i++) {
        if ((pBP[i] < 0) && (pBP[i] > -4)) {
            bonus_cnt++;
            if ((pBP[i] == -3) && (pStructure[i - 1] == ')')) {
                bonus++;
            }
            if ((pBP[i] == -2) && (pStructure[i - 1] == '(')) {
                bonus++;
            }
            if ((pBP[i] == -1) && (pStructure[i - 1] != '.')) {
                bonus++;
            }
        }

        if (pBP[i] > i) {
            int l;
            bonus_cnt++;
            for (l = 1; l <= pOpts.mBasePair[0].i; l++) {
                if ((i == pOpts.mBasePair[l].i) && (mBP[i] == pOpts.mBasePair[l].j)) {
                    bonus++;
                }
            }
        }
    }

    if (bonus_cnt > bonus) {
        std::cerr << "Error: Could not enforce all constraints." << std::endl;
    }
    bonus *= BONUS;

    pEnergy += bonus; /*remove bonus energies from result */

    if (pOpts.mBacktrackType == 'C') {
        return (double) mArrayC.mC[mIdex[pStrLen] + 1] / 100.0;
    } else if (pOpts.mBacktrackType == 'M') {
        return (double) mArrayML.mFML[mIdex[pStrLen] + 1] / 100.0;
    } else {
        return (double) pEnergy / 100.0;
    }
}

//
// VR16FoldArrayC methods
//
int VR16FoldArrayC::init_arrays(int pStrLen) {
    if (pStrLen < 1) {
        std::cerr << "Error: Length must be greater 0" << std::endl;
        return 1;
    }

    if (pStrLen > mInitLength) {
        clear_arrays();
        resize_arrays(pStrLen);
    }

    for (int j = 1; j <= pStrLen; ++j) {
        for (int i = (j > TURN ? (j - TURN) : 1); i < j; ++i) {
            mC[mIdex[j] + i] = mParams.INF;
        }
    }

    for (int j = 1; j <= pStrLen; ++j) {
        mDMLi[j] = mParams.INF;
        mDMLi1[j] = mParams.INF;
        mDMLi2[j] = mParams.INF;
    }

    mIdDMLi2 = 2;
    mIdDMLi1 = 1;
    mIdDMLi = 0;

    mIdCC1 = 1;
    mIdCC = 0;

    mInitLength = pStrLen;

    return 0;
}

void VR16FoldArrayC::clear_arrays() {
    mC.clear();
    mCC.clear();
    mCC1.clear();

    mDMLi.clear();
    mDMLi1.clear();
    mDMLi2.clear();
}

void VR16FoldArrayC::resize_arrays(int pStrLen) {
    mC.resize((unsigned) (pStrLen * (pStrLen + 1)) / 2 + 2);
    mCC.resize((unsigned) pStrLen + 2);
    mCC1.resize((unsigned) pStrLen + 2);

    mDMLi.resize((unsigned) pStrLen + 1);
    mDMLi1.resize((unsigned) pStrLen + 1);
    mDMLi2.resize((unsigned) pStrLen + 1);
}

int VR16FoldArrayC::get_init_bonus(int pI, int pJ, std::vector<int> &pBP) {
    int bonus = 0;

    /* enforcing structure constraints */
    if ((pBP[pI] == pJ) || (pBP[pI] == -1) || (pBP[pI] == -2)) {
        bonus -= BONUS;
    }
    if ((pBP[pJ] == -1) || (pBP[pJ] == -3)) {
        bonus -= BONUS;
    }

    return bonus;
}

void VR16FoldArrayC::fill_arrays(
        int pI,
        int pJ,
        std::vector<int> &pBP,
        std::vector<char> &pPtype,
        std::vector<int> &pS1,
        std::vector<int> &pFML,
        std::string &pString) {
    int strLen = (int) pString.size();
    int ij = mIdex[pJ] + pI;
    int bonus = get_init_bonus(pI, pJ, pBP);
    int type_1 = get_type_1(pI, pJ, strLen, pPtype, pBP);
    int no_close = (((type_1 == 3) || (type_1 == 4)) && mOpts.mNoClosingGU && (bonus == 0));

    if (type_1) { /* we have a pair */
        int new_c = 0;
        int stackEnergy = mParams.INF;

        /* hairpin ----------------------------------------------*/
        if (no_close) {
            new_c = FORBIDDEN;
        } else {
            new_c = hairpin_e(pJ - pI - 1, type_1, pS1[pI + 1], pS1[pJ - 1], pI - 1, pString);
        }

        /* multi-loop decomposition ------------------------*/
        if (!no_close) {
            int ml_energy = calc_ml_decomp(pI, pJ, type_1, pS1);
            if (ml_energy < new_c) {
                new_c = ml_energy;
            }
        }

        /* coaxial stacking of (i.j) with (i+1.k) or (k+1.j-1) */
        if (mOpts.mDangles == 3) {
            int decomp = calc_coaxial_stack(pI, pJ, type_1, pPtype, pFML);
            if (decomp < new_c) {
                new_c = decomp;
            }
        }

        /*--------------------------------------------------------
         check for elementary structures involving more than one
         closing pair.
         --------------------------------------------------------*/
        calc_elementary_struct(pI, pJ, type_1, no_close, pPtype, pS1, new_c, stackEnergy);

        set_c(pI, pJ, bonus, new_c, stackEnergy);

    } /* end >> if (pair) << */
    else {
        mC[ij] = mParams.INF;
    }

}

int VR16FoldArrayC::get_type_1(int pI, int pJ, int pStrLen, std::vector<char> &pPtype, std::vector<int> &pBP) {
    int type_1 = pPtype[mIdex[pJ] + pI];
    int max_separation = (int) ((1.0 - LOCALITY) * (double) (pStrLen - 2)); /* not in use */

    if ((pBP[pI] == -4) || (pBP[pJ] == -4)) {
        type_1 = 0;
    }

    if (pJ - pI - 1 > max_separation) {
        type_1 = 0; /* forces locality degree */
    }

    return type_1;
}

int VR16FoldArrayC::calc_ml_decomp(int pI, int pJ, int pType1, std::vector<int> &pS1) {
    int ml_energy;
    int decomp = get_dmlx_val(mIdDMLi1, pJ - 1);
    int decomp_t;

    if (mOpts.mDangles) {
        int d3 = 0;
        int d5 = 0;
        int tt = mPairMat.mRtype[pType1];
        d3 = mParams.mDangle3[tt][pS1[pI + 1]];
        d5 = mParams.mDangle5[tt][pS1[pJ - 1]];
        if (mOpts.mDangles == 2) /* double dangles */
        {
            decomp += d5 + d3;
        } else { /* normal dangles */
            decomp_t = get_dmlx_val(mIdDMLi2, pJ - 1) + d3 + mParams.mMLBase;
            if (decomp_t < decomp) {
                decomp = decomp_t;
            }

            decomp_t = get_dmlx_val(mIdDMLi1, pJ - 2) + d5 + mParams.mMLBase;
            if (decomp_t < decomp) {
                decomp = decomp_t;
            }

            decomp_t = get_dmlx_val(mIdDMLi2, pJ - 2) + d5 + d3 + 2 * mParams.mMLBase;
            if (decomp_t < decomp) {
                decomp = decomp_t;
            }
        }
    }

    ml_energy = mParams.mMLclosing + mParams.mMLintern[pType1] + decomp;

    return ml_energy;
}

int VR16FoldArrayC::calc_coaxial_stack(
        int pI,
        int pJ,
        int pType1,
        std::vector<char> &pPtype,
        std::vector<int> &pFML) {
    int type_2;
    int decomp = mParams.INF;
    int decomp_t;

    for (int k = pI + 2 + TURN; k < pJ - 2 - TURN; ++k) {
        type_2 = pPtype[mIdex[k] + pI + 1];
        type_2 = mPairMat.mRtype[type_2];
        if (type_2) {
            decomp_t = mC[mIdex[k] + pI + 1] + mParams.mStack[pType1][type_2] + pFML[mIdex[pJ - 1] + k + 1];
            if (decomp_t < decomp) {
                decomp = decomp_t;
            }
        }
        type_2 = pPtype[mIdex[pJ - 1] + k + 1];
        type_2 = mPairMat.mRtype[type_2];
        if (type_2) {
            decomp_t = mC[mIdex[pJ - 1] + k + 1] + mParams.mStack[pType1][type_2] + pFML[mIdex[k] + pI + 1];
            if (decomp_t < decomp) {
                decomp = decomp_t;
            }
        }
    }
    /* no TermAU penalty if coax stack */
    decomp += 2 * mParams.mMLintern[1] + mParams.mMLclosing;

    return decomp;
}

void VR16FoldArrayC::calc_elementary_struct(
        int pI,
        int pJ,
        int pType1,
        int pNoClose,
        std::vector<char> &pPtype,
        std::vector<int> &pS1,
        int &pNewC,
        int &pStackEnergy) {
    int type_2;
    int energy;
    int minq;

    int max_p = pJ - 2 - TURN;
    if (pI + mParams.MAXLOOP + 1 < max_p) {
        max_p = pI + mParams.MAXLOOP + 1;
    }
    for (int p = max_p; p >= pI + 1; --p) {
        minq = pJ - pI + p - mParams.MAXLOOP - 2;
        if (minq < p + 1 + TURN) {
            minq = p + 1 + TURN;
        }
        for (int q = minq; q < pJ; ++q) {
            type_2 = pPtype[mIdex[q] + p];
            if (type_2 == 0) {
                continue;
            }

            if (pNewC < mC[mIdex[q] + p] + mPpIL.mMinIL[p - pI - 1][pJ - q - 1]) {
                continue;
            }

            type_2 = mPairMat.mRtype[type_2];
            if (mOpts.mNoClosingGU) {
                if (pNoClose || (type_2 == 3) || (type_2 == 4)) {
                    if ((p > pI + 1) || (q < pJ - 1)) {
                        continue; /* continue unless stack */
                    }
                }
            }

            energy = mPpIL.loop_energy(p - pI - 1, pJ - q - 1, pType1, type_2, pS1[pI + 1], pS1[pJ - 1],
                                       pS1[p - 1], pS1[q + 1]);
            if (energy + mC[mIdex[q] + p] < pNewC) {
                pNewC = energy + mC[mIdex[q] + p];
            }

            if ((p == pI + 1) && (pJ == q + 1)) {
                pStackEnergy = energy; /* remember stack energy */
            }

        } /* end q-loop */
    } /* end p-loop */
}

void VR16FoldArrayC::set_c(
        int pI,
        int pJ,
        int pBonus,
        int pNewC,
        int pStackEnergy) {
    int energy_t = get_ccx_val(mIdCC1, pJ - 1) + pStackEnergy;
    int ij = mIdex[pJ] + pI;

    if (energy_t < pNewC) {
        pNewC = energy_t;
    }

    set_ccx_val(mIdCC, pJ, pNewC + pBonus);
    if (mOpts.mNoLonelyPairs) {
        mC[ij] = get_ccx_val(mIdCC1, pJ - 1) + pStackEnergy + pBonus;
    } else {
        mC[ij] = get_ccx_val(mIdCC, pJ);
    }
}

int VR16FoldArrayC::hairpin_e(int pSize, int pType, int pSi1, int pSj1, int pIdx, std::string &pString) {
    int energy;
    std::string subStr6;
    std::string subStr5;

    energy = (pSize <= 30) ? mParams.mHairpin[pSize] :
             mParams.mHairpin[30] + (int) (mParams.mlxc * std::log((pSize) / 30.0));

    if (mOpts.mTetraLoop) {
        if (pSize == 4) { /* check for tetraloop bonus */
            subStr6 = pString.substr((unsigned) pIdx, 6);
            std::map<std::string, int>::const_iterator found = mEn.mTetraloops.find(subStr6);
            if (found != mEn.mTetraloops.end()) {
                energy += mParams.mTetraEnergy[(*found).second];
            }
        }
    }

    if (pSize == 3) {
        if (mOpts.mTriLoop) {
            subStr5 = pString.substr((unsigned) pIdx, 5);
            std::map<std::string, int>::const_iterator found = mEn.mTriloops.find(subStr5);
            if (found != mEn.mTriloops.end()) {
                energy += mParams.mTriloopE[(*found).second];
            }
        }

        if (pType > 2) /* neither CG nor GC */
        {
            energy += mParams.mTerminalAU; /* penalty for closing AU GU pair */
        }
    } else {
        /* no mismatches for tri-loops */
        energy += mParams.mMismatchH[pType][pSi1][pSj1];
    }

    return energy;
}

void VR16FoldArrayC::rotate_aux_idexes(int pStrLen) {
    /* rotate the auxilliary arrays */
    int idTmp;

    idTmp = mIdDMLi2;
    mIdDMLi2 = mIdDMLi1;
    mIdDMLi1 = mIdDMLi;
    mIdDMLi = idTmp;

    idTmp = mIdCC1;
    mIdCC1 = mIdCC;
    mIdCC = idTmp;

    for (int j = 1; j <= pStrLen; ++j) {
        set_ccx_val(mIdCC, j, mParams.INF);
        set_dmlx_val(mIdDMLi, j, mParams.INF);
    }
}

//
// VR16FoldArrayML methods
//
int VR16FoldArrayML::init_arrays(int pStrLen) {
    if (pStrLen < 1) {
        std::cerr << "Error: Length must be greater 0" << std::endl;
        return 1;
    }

    if (pStrLen > mInitLength) {
        clear_arrays();
        resize_arrays(pStrLen);
    }

    for (int j = 1; j <= pStrLen; ++j) {
        mFmi[j] = mParams.INF;
    }

    for (int j = 1; j <= pStrLen; ++j) {
        for (int i = (j > TURN ? (j - TURN) : 1); i < j; ++i) {
            mFML[mIdex[j] + i] = mParams.INF;
            if (mUniqML) {
                mFM1[mIdex[j] + i] = mParams.INF;
            }
        }
    }

    mInitLength = pStrLen;

    return 0;
}

void VR16FoldArrayML::clear_arrays() {
    mFML.clear();
    mFM1.clear();
    mFmi.clear();
    mInitLength = 0;
}

void VR16FoldArrayML::resize_arrays(int pStrLen) {
    mFML.resize((unsigned) (pStrLen * (pStrLen + 1)) / 2 + 2);
    if (mUniqML) {
        mFM1.resize((unsigned) (pStrLen * (pStrLen + 1)) / 2 + 2);
    }
    mFmi.resize((unsigned) pStrLen + 1);

}

void VR16FoldArrayML::fill_arrays(
        int pI,
        int pJ,
        int pStrLen,
        std::vector<char> &pPtype,
        std::vector<int> &pS1,
        VR16FoldArrayC &pArrayC) {
    int ij = mIdex[pJ] + pI;
    int type_1 = pPtype[ij];

    int new_fML;
    int fML_t;
    int energy;
    int decomp;

    new_fML = mFML[ij + 1] + mParams.mMLBase;
    fML_t = mFML[mIdex[pJ - 1] + pI] + mParams.mMLBase;
    if (fML_t < new_fML) {
        new_fML = fML_t;
    }

    energy = pArrayC.mC[ij] + mParams.mMLintern[type_1];
    if (mOpts.mDangles == 2) { /* double dangles */
        energy += (pI == 1) ? /* works also for circfold */
                  mParams.mDangle5[type_1][pS1[pStrLen]] : mParams.mDangle5[type_1][pS1[pI - 1]];
        /* if (j<length) */
        energy += mParams.mDangle3[type_1][pS1[pJ + 1]];
    }
    if (energy < new_fML) {
        new_fML = energy;
    }

    if (mUniqML) {
        int energy_t = mFM1[mIdex[pJ - 1] + pI] + mParams.mMLBase;
        if (energy_t < energy) {
            mFM1[ij] = energy_t;
        }
    }

    if (mOpts.mDangles % 2 == 1) { /* normal dangles */
        int tt = pPtype[ij + 1]; /* i+1,j */
        fML_t = pArrayC.mC[ij + 1] + mParams.mDangle5[tt][pS1[pI]] + mParams.mMLintern[tt];
        fML_t += mParams.mMLBase;
        if (fML_t < new_fML) {
            new_fML = fML_t;
        }

        tt = pPtype[mIdex[pJ - 1] + pI];
        fML_t = pArrayC.mC[mIdex[pJ - 1] + pI] + mParams.mDangle3[tt][pS1[pJ]] + mParams.mMLintern[tt];
        fML_t += mParams.mMLBase;
        if (fML_t < new_fML) {
            new_fML = fML_t;
        }

        tt = pPtype[mIdex[pJ - 1] + pI + 1];
        fML_t = pArrayC.mC[mIdex[pJ - 1] + pI + 1] + mParams.mDangle5[tt][pS1[pI]]
                + mParams.mDangle3[tt][pS1[pJ]] + mParams.mMLintern[tt] + 2 * mParams.mMLBase;
        if (fML_t < new_fML) {
            new_fML = fML_t;
        }
    }

    /* modular decomposition -------------------------------*/
    decomp = mParams.INF;
    for (int k = pI + 1 + TURN; k <= pJ - 2 - TURN; ++k) {
        int decomp_t = mFmi[k] + mFML[mIdex[pJ] + k + 1];
        if (decomp_t < decomp) {
            decomp = decomp_t;
        }
    }
    pArrayC.set_dmli_val(pJ, decomp); /* store for use in ML decompositon */
    if (decomp < new_fML) {
        new_fML = decomp;
    }

    /* coaxial stacking */
    if (mOpts.mDangles == 3) {
        /* additional ML decomposition as two coaxially stacked helices */
        decomp = mParams.INF;
        for (int k = pI + 1 + TURN; k <= pJ - 2 - TURN; ++k) {
            type_1 = pPtype[mIdex[k] + pI];
            type_1 = mPairMat.mRtype[type_1];
            int type_2 = pPtype[mIdex[pJ] + k + 1];
            type_2 = mPairMat.mRtype[type_2];
            if (type_1 && type_2) {
                int decomp_t = pArrayC.mC[mIdex[k] + pI] + pArrayC.mC[mIdex[pJ] + k + 1]
                               + mParams.mStack[type_1][type_2];
                if (decomp_t < decomp) {
                    decomp = decomp_t;
                }
            }
        }

        decomp += 2 * mParams.mMLintern[1]; /* no TermAU penalty if coax stack */
        if (decomp < new_fML) {
            new_fML = decomp;
        }
    }

    mFML[ij] = mFmi[pJ] = new_fML; /* substring energy */
}

void VR16FoldArrayML::reset_fmi(int pStrLen) {
    for (int j = 1; j <= pStrLen; ++j) {
        mFmi[j] = mParams.INF;
    }
}

//
// VR16FoldBackTrack methods
//
void VR16FoldBackTrack::backtrack(
        std::string &pString,
        int pS,
        std::vector<char> &pPtype,
        std::vector<int> &pS1,
        std::vector<int> &pBP,
        VR16FoldArrayC &pArrayC,
        VR16FoldArrayML &pArrayML,
        std::vector<int> &pF5) {

    /*------------------------------------------------------------------
     trace back through the "c", "f5" and "fML" arrays to get the
     base pairing list. No search for equivalent structures is done.
     This is fast, since only few structure elements are recalculated.
     ------------------------------------------------------------------*/

    /* normally s=0.
     If s>0 then s items have been already pushed onto the sector stack */
//    int i, j, k, length, energy, new;
    int no_close, type_1, type_2, tt;
    int bonus;
    int b = 0;
    int strLen;
    int i;
    int j;
    int k;
    int energy;
    int newE;

    strLen = (int) pString.size();

    if (pS == 0) {
        mSector[++pS].i = 1;
        mSector[pS].j = strLen;
        mSector[pS].ml = (mOpts.mBacktrackType == 'M') ? 1 : ((mOpts.mBacktrackType == 'C') ? 2 : 0);
    }

    while (pS > 0) {
        int ml, fij, fi, traced, i1, j1, d3, d5, mm, p, q, jj = 0;
        int cij = 0;
        int canonical = 1; /* (i,j) closes a canonical structure */
        i = mSector[pS].i;
        j = mSector[pS].j;
        ml = mSector[pS--].ml; /* ml is a flag indicating if backtracking is to
         occur in the fML- (1) or in the f-array (0) */
        if (ml == 2) {
            mOpts.mBasePair[++b].i = i;
            mOpts.mBasePair[b].j = j;
            goto repeat1;
        }

        if (j < i + TURN + 1) {
            continue; /* no more pairs in this interval */
        }

        fij = (ml) ? pArrayML.mFML[mIdex[j] + i] : pF5[j];
        fi = (ml) ? (pArrayML.mFML[mIdex[j - 1] + i] + mParams.mMLBase) : pF5[j - 1];

        if (fij == fi) { /* 3' end is unpaired */
            mSector[++pS].i = i;
            mSector[pS].j = j - 1;
            mSector[pS].ml = ml;
            continue;
        }

        if (ml == 0) { /* backtrack in f5 */
            /* j or j-1 is paired. Find pairing partner */
            for (k = j - TURN - 1, traced = 0; k >= 1; k--) {
                int cc, en;
                jj = k - 1;
                type_1 = pPtype[mIdex[j - 1] + k];
                if ((type_1) && (mOpts.mDangles % 2 == 1)) {
                    cc = pArrayC.mC[mIdex[j - 1] + k] + mParams.mDangle3[type_1][pS1[j]];
                    if (type_1 > 2) {
                        cc += mParams.mTerminalAU;
                    }
                    if (fij == cc + pF5[k - 1]) {
                        traced = j - 1;
                    }
                    if (k > i) {
                        if (fij == pF5[k - 2] + cc + mParams.mDangle5[type_1][pS1[k - 1]]) {
                            traced = j - 1;
                            jj = k - 2;
                        }
                    }
                }
                type_1 = pPtype[mIdex[j] + k];
                if (type_1) {
                    cc = pArrayC.mC[mIdex[j] + k];
                    if (type_1 > 2) {
                        cc += mParams.mTerminalAU;
                    }
                    en = cc + pF5[k - 1];
                    if (mOpts.mDangles == 2) {
                        if (k > 1) {
                            en += mParams.mDangle5[type_1][pS1[k - 1]];
                        }
                        if (j < strLen) {
                            en += mParams.mDangle3[type_1][pS1[j + 1]];
                        }
                    }
                    if (fij == en) {
                        traced = j;
                    }
                    if ((mOpts.mDangles % 2 == 1) && (k > 1)) {
                        if (fij == pF5[k - 2] + cc + mParams.mDangle5[type_1][pS1[k - 1]]) {
                            traced = j;
                            jj = k - 2;
                        }
                    }
                }
                if (traced) {
                    break;
                }
            }

            if (!traced) {
                std::cerr << "Error: Backtrack failed in f5s." << std::endl;
            }
            mSector[++pS].i = 1;
            mSector[pS].j = jj;
            mSector[pS].ml = ml;

            i = k;
            j = traced;
            mOpts.mBasePair[++b].i = i;
            mOpts.mBasePair[b].j = j;
            goto repeat1;
        } else { /* trace back in fML array */
            int cij1 = mParams.INF;
            int ci1j = mParams.INF;
            int ci1j1 = mParams.INF;
            if (pArrayML.mFML[mIdex[j] + i + 1] + mParams.mMLBase == fij) { /* 5' end is unpaired */
                mSector[++pS].i = i + 1;
                mSector[pS].j = j;
                mSector[pS].ml = ml;
                continue;
            }

            tt = pPtype[mIdex[j] + i];
            cij = pArrayC.mC[mIdex[j] + i] + mParams.mMLintern[tt];
            if (mOpts.mDangles == 2) { /* double dangles, works also for circfold */
                cij += (i == 1) ? mParams.mDangle5[tt][pS1[strLen]] : mParams.mDangle5[tt][pS1[i - 1]];
                /* if (j<length) */
                cij += mParams.mDangle3[tt][pS1[j + 1]];
            } else if (mOpts.mDangles % 2 == 1) { /* normal dangles */
                tt = pPtype[mIdex[j] + i + 1];
                ci1j = pArrayC.mC[mIdex[j] + i + 1] + mParams.mDangle5[tt][pS1[i]] + mParams.mMLintern[tt] +
                       mParams.mMLBase;

                tt = pPtype[mIdex[j - 1] + i];
                cij1 = pArrayC.mC[mIdex[j - 1] + i] + mParams.mDangle3[tt][pS1[j]] + mParams.mMLintern[tt] +
                       mParams.mMLBase;

                tt = pPtype[mIdex[j - 1] + i + 1];
                ci1j1 = pArrayC.mC[mIdex[j - 1] + i + 1] + mParams.mDangle5[tt][pS1[i]] +
                        mParams.mDangle3[tt][pS1[j]]
                        + mParams.mMLintern[tt] + 2 * mParams.mMLBase;
            }

            if ((fij == cij) || (fij == ci1j) || (fij == cij1) || (fij == ci1j1)) {
                /* found a pair */
                if (fij == ci1j) {
                    i++;
                } else if (fij == cij1) {
                    j--;
                } else if (fij == ci1j1) {
                    i++;
                    j--;
                }
                mOpts.mBasePair[++b].i = i;
                mOpts.mBasePair[b].j = j;
                goto repeat1;
            }

            for (k = i + 1 + TURN; k <= j - 2 - TURN; k++) {
                if (fij == (pArrayML.mFML[mIdex[k] + i] + pArrayML.mFML[mIdex[j] + k + 1])) {
                    break;
                }
            }

            if ((mOpts.mDangles == 3) && (k > j - 2 - TURN)) { /* must be coax stack */
                ml = 2;
                for (k = i + 1 + TURN; k <= j - 2 - TURN; k++) {
                    type_1 = pPtype[mIdex[k] + i];
                    type_1 = mPairMat.mRtype[type_1];
                    type_2 = pPtype[mIdex[j] + k + 1];
                    type_2 = mPairMat.mRtype[type_2];
                    if (type_1 && type_2) {
                        if (fij == pArrayC.mC[mIdex[k] + i] + pArrayC.mC[mIdex[j] + k + 1]
                                   + mParams.mStack[type_1][type_2] + 2 * mParams.mMLintern[1]) {
                            break;
                        }
                    }
                }
            }

            mSector[++pS].i = i;
            mSector[pS].j = k;
            mSector[pS].ml = ml;
            mSector[++pS].i = k + 1;
            mSector[pS].j = j;
            mSector[pS].ml = ml;

            if (k > j - 2 - TURN) {
                std::cerr << "Error: backtrack failed in fML." << std::endl;
            }

            continue;
        }

        repeat1:

        /*----- begin of "repeat:" -----*/
        if (canonical) {
            cij = pArrayC.mC[mIdex[j] + i];
        }

        type_1 = pPtype[mIdex[j] + i];

        bonus = 0;

        if ((pBP[i] == j) || (pBP[i] == -1) || (pBP[i] == -2)) {
            bonus -= BONUS;
        }
        if ((pBP[j] == -1) || (pBP[j] == -3)) {
            bonus -= BONUS;
        }

        if (mOpts.mNoLonelyPairs) {
            if (cij == pArrayC.mC[mIdex[j] + i]) {
                /* (i.j) closes canonical structures, thus
                 (i+1.j-1) must be a pair                */
                type_2 = pPtype[mIdex[j - 1] + i + 1];
                type_2 = mPairMat.mRtype[type_2];
                cij -= mParams.mStack[type_1][type_2] + bonus;
                mOpts.mBasePair[++b].i = i + 1;
                mOpts.mBasePair[b].j = j - 1;
                i++;
                j--;
                canonical = 0;
                goto repeat1;
            }
        }

        canonical = 1;

        no_close = (((type_1 == 3) || (type_1 == 4)) && mOpts.mNoClosingGU && (bonus == 0));
        if (no_close) {
            if (cij == FORBIDDEN) {
                continue;
            }
        } else if (cij == pArrayC.hairpin_e(j - i - 1, type_1, pS1[i + 1], pS1[j - 1], i - 1, pString) + bonus) {
            continue;
        }

        int max_p = (((j - 2 - TURN) < (i + mParams.MAXLOOP + 1)) ? (j - 2 - TURN) : (i + mParams.MAXLOOP + 1));
        for (p = i + 1; p <= max_p; p++) {
            int minq;
            minq = j - i + p - mParams.MAXLOOP - 2;
            if (minq < p + 1 + TURN) {
                minq = p + 1 + TURN;
            }
            for (q = j - 1; q >= minq; q--) {

                type_2 = pPtype[mIdex[q] + p];
                if (type_2 == 0) {
                    continue;
                }
                type_2 = mPairMat.mRtype[type_2];
                if (mOpts.mNoClosingGU) {
                    if (no_close || (type_2 == 3) || (type_2 == 4)) {
                        if ((p > i + 1) || (q < j - 1)) {
                            continue; /* continue unless stack */
                        }
                    }
                }

                energy = mPpIL.loop_energy(p - i - 1, j - q - 1, type_1, type_2, pS1[i + 1], pS1[j - 1],
                                           pS1[p - 1], pS1[q + 1]);

                newE = energy + pArrayC.mC[mIdex[q] + p] + bonus;
                traced = (cij == newE);
                if (traced) {
                    mOpts.mBasePair[++b].i = p;
                    mOpts.mBasePair[b].j = q;
                    i = p, j = q;
                    goto repeat1;
                }
            }
        }

        /* end of repeat: --------------------------------------------------*/

        /* (i.j) must close a multi-loop */
        tt = mPairMat.mRtype[type_1];
        mm = bonus + mParams.mMLclosing + mParams.mMLintern[tt];
        d5 = mParams.mDangle5[tt][pS1[j - 1]];
        d3 = mParams.mDangle3[tt][pS1[i + 1]];
        i1 = i + 1;
        j1 = j - 1;
        mSector[pS + 1].ml = 1;
        mSector[pS + 2].ml = 1;

        for (k = i + 2 + TURN; k < j - 2 - TURN; k++) {
            int en;
            en = pArrayML.mFML[mIdex[k] + i + 1] + pArrayML.mFML[mIdex[j - 1] + k + 1] + mm;
            if (mOpts.mDangles == 2) /* double dangles */
            {
                en += d5 + d3;
            }
            if (cij == en) {
                break;
            }
            if (mOpts.mDangles % 2 == 1) { /* normal dangles */
                if (cij == (pArrayML.mFML[mIdex[k] + i + 2] + pArrayML.mFML[mIdex[j - 1] + k + 1] + mm + d3 +
                            mParams.mMLBase)) {
                    i1 = i + 2;
                    break;
                }
                if (cij == (pArrayML.mFML[mIdex[k] + i + 1] + pArrayML.mFML[mIdex[j - 2] + k + 1] + mm + d5 +
                            mParams.mMLBase)) {
                    j1 = j - 2;
                    break;
                }
                if (cij == (pArrayML.mFML[mIdex[k] + i + 2] + pArrayML.mFML[mIdex[j - 2] + k + 1]
                            + mm + d3 + d5 + mParams.mMLBase + mParams.mMLBase)) {
                    i1 = i + 2;
                    j1 = j - 2;
                    break;
                }
            }
            /* coaxial stacking of (i.j) with (i+1.k) or (k.j-1) */
            /* use MLintern[1] since coax stacked pairs don't get TerminalAU */
            if (mOpts.mDangles == 3) {
                type_2 = pPtype[mIdex[k] + i + 1];
                type_2 = mPairMat.mRtype[type_2];
                if (type_2) {
                    en = pArrayC.mC[mIdex[k] + i + 1] + mParams.mStack[type_1][type_2] +
                         pArrayML.mFML[mIdex[j - 1] + k + 1];
                    if (cij == en + 2 * mParams.mMLintern[1] + mParams.mMLclosing) {
                        ml = 2;
                        mSector[pS + 1].ml = 2;
                        break;
                    }
                }
                type_2 = pPtype[mIdex[j - 1] + k + 1];
                type_2 = mPairMat.mRtype[type_2];
                if (type_2) {
                    en = pArrayC.mC[mIdex[j - 1] + k + 1] + mParams.mStack[type_1][type_2] +
                         pArrayML.mFML[mIdex[k] + i + 1];
                    if (cij == en + 2 * mParams.mMLintern[1] + mParams.mMLclosing) {
                        mSector[pS + 2].ml = 2;
                        break;
                    }
                }
            }

        }
        if (k <= j - 3 - TURN) { /* found the decomposition */
            mSector[++pS].i = i1;
            mSector[pS].j = k;
            mSector[++pS].i = k + 1;
            mSector[pS].j = j1;
        } else {
            std::cerr << "Error: backtracking failed in repeat." << std::endl;
        }

    }

    mOpts.mBasePair[0].i = b; /* save the total number of base pairs */
}

void VR16FoldBackTrack::parenthesis_structure(std::string &pStructure, int pStrLen) {
    int n, k;

    for (n = 0; n <= pStrLen - 1; pStructure[n++] = '.');
    pStructure[pStrLen] = '\0';

    for (k = 1; k <= mOpts.mBasePair[0].i; k++) {
        pStructure[mOpts.mBasePair[k].i - 1] = '(';
        pStructure[mOpts.mBasePair[k].j - 1] = ')';
    }
}

} // namespace vr16
