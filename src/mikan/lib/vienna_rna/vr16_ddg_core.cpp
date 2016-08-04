#include <mikan/lib/vienna_rna/include/vr16_ddg_core.hpp>
#include <iostream>

namespace vr16
{

//
// VR16DuplexRet methods
//
void VR16DuplexRet::init_ret_vals(int pSize)
{
    if (pSize == mVecSize)
    {
        return;
    }

    mStructure.resize((unsigned)pSize);
    mL1.resize((unsigned)pSize);
    mStartMiR.resize((unsigned)pSize);
    mEndMiR.resize((unsigned)pSize);
    mStartTarget.resize((unsigned)pSize);
    mEndTarget.resize((unsigned)pSize);
    mMiRLen.resize((unsigned)pSize);
    mTargetLen.resize((unsigned)pSize);
    mDGall.resize((unsigned)pSize);
    mDG5.resize((unsigned)pSize);
    mDG3.resize((unsigned)pSize);

    mVecSize = pSize;
}

//
// VR16DDG4Ret methods
//
void VR16DDG4Ret::init_ret_vals(int pSize)
{
    if (pSize == mVecSize)
    {
        return;
    }

    mDG0.resize((unsigned)pSize);
    mD1Array.resize((unsigned)pSize);
    mD2Array.resize((unsigned)pSize);
    mDDGArray.resize((unsigned)pSize);
    for (int i = 0; i < pSize; ++i)
    {
        mD1Array[i].resize(ARRSIZE);
        mD2Array[i].resize(ARRSIZE);
        mDDGArray[i].resize(ARRSIZE);
    }
    mDDGSum.resize((unsigned)pSize);
    mP.resize((unsigned)pSize);

    mVecSize = pSize;
}

//
// VR16DDGWorkSpace methods
//
void VR16DDGWorkSpace::init_workspace()
{
    mPairMat.make_pair_matrix(mDupOpts.mEnergySet, mDupOpts.mNonStandards, mDupOpts.mNoGU);
}

void VR16DDGWorkSpace::preppare_duplexfold(int pSize)
{
    mDupRet.init_ret_vals(pSize);
}

int VR16DDGWorkSpace::duplexfold(
        int pRetIdx,
        std::string &pS1,
        std::string &pS2,
        std::vector<int> &pArrayI,
        std::vector<int> &pArrayJ)
{
    int retVal;
    int l2;

    retVal = mDup.duplexfold(pS1, pS2, pArrayI, pArrayJ);
    if (retVal != 0)
    {
        return retVal;
    }

    l2 =  (int)mDup.mStructure.size() - mDup.mL1 - 1;

    mDupRet.mStructure[pRetIdx] = mDup.mStructure;
    mDupRet.mL1[pRetIdx] = mDup.mL1;
    mDupRet.mStartMiR[pRetIdx] = mDup.mI + 1 - mDup.mL1;
    mDupRet.mEndMiR[pRetIdx] = mDup.mI;
    mDupRet.mStartTarget[pRetIdx] = mDup.mJ;
    mDupRet.mEndTarget[pRetIdx] = mDup.mJ + l2 - 1;
    mDupRet.mMiRLen[pRetIdx] = mDup.mL1;
    mDupRet.mTargetLen[pRetIdx] = l2;
    mDupRet.mDGall[pRetIdx] = mDup.mEnergy;
    mDupRet.mDG5[pRetIdx] = mDup.mE5p / 100.0;
    mDupRet.mDG3[pRetIdx] = mDup.mE3p / 100.0;

    return 0;
}

void VR16DDGWorkSpace::print_duplexfold_ret_vals(int pRetIdx)
{
    std::cout << "Structure:     " << mDupRet.mStructure[pRetIdx] << std::endl;
    std::cout << "Start miRNA:   " << mDupRet.mStartMiR[pRetIdx] << std::endl;
    std::cout << "End miRNA:     " << mDupRet.mEndMiR[pRetIdx] << std::endl;
    std::cout << "Start target:  " << mDupRet.mStartTarget[pRetIdx] << std::endl;
    std::cout << "End target:    " << mDupRet.mEndTarget[pRetIdx] << std::endl;
    std::cout << "miRNA length:  " << mDupRet.mMiRLen[pRetIdx] << std::endl;
    std::cout << "Target length: " << mDupRet.mTargetLen[pRetIdx] << std::endl;
    std::cout << "dGall:         " << mDupRet.mDGall[pRetIdx] << std::endl;
    std::cout << "dG5:           " << mDupRet.mDG5[pRetIdx] << std::endl;
    std::cout << "dG3:           " << mDupRet.mDG3[pRetIdx] << std::endl;
}

void VR16DDGWorkSpace::prepare_ddg4(
        int pSize,
        int pUpRest,
        int pTargetStart,
        int pRestrictedFrom,
        int pRestrictedTo)
{
    mDDG4Ret.init_ret_vals(pSize);

    mUpRest = pUpRest;
    mTargetStart = pTargetStart;
    mRestrictedFrom = pRestrictedFrom;
    mRestrictedTo = pRestrictedTo;

    mDG4Opts.mDoBacktrack = false;
    mDG4Opts.mDangles = 2;
    mDG4Opts.mFoldConstrained = false;
}

int VR16DDGWorkSpace::calc_ddg4(int pRetIdx, std::string &pString)
{
    double min_en0;
    double min_en1;
    std::string struct0;
    std::string struct1;
    double dg0;
    double dg1;
    double dg2 = 0.0;
    double kT;
    double ddGsum;
    std::vector<double> ddG((unsigned)mDDG4Ret.ARRSIZE);
    std::string constrains0(pString.size(), 0);

    /* Fold with no constraint ********************************************/
    mDG4Opts.mFoldConstrained = false;
    struct0 = constrains0;
    min_en0 = mFold.fold(pString, struct0);

    kT = (mTemperature + 273.15) * 1.98717 / 1000.;     /* in Kcal */
    mDG4Opts.mPfScale = std::exp(-1.0 * (SFACT * min_en0) / kT / pString.size());
    struct0 = constrains0;
    dg0 = mPf.pf_fold(pString, struct0, mTemperature);
    mDDG4Ret.mDG0[pRetIdx] = dg0;

    mDG4Opts.mFoldConstrained = true;
    mDDG4Ret.mNExps = 0;
    for (int i = mRestrictedFrom; i <= mRestrictedTo; ++i)
    {
        /* Build constraint string */
        std::string constrains1(pString.size(), ' ');
        for (int j = (mTargetStart + mUpRest); j > (mTargetStart - i); --j)
        {
            constrains1[j] = 'x';
        }

        if (i < mRestrictedTo)
        {
            constrains1[mTargetStart - i] = '|';
        }

        /* Fold with constraint *******************************************/
        struct1 = constrains1;
        min_en1 = mFold.fold(pString, struct1);

        mDG4Opts.mPfScale = std::exp(-1.0 * (SFACT * min_en1) / kT / pString.size());
        struct1 = constrains1;
        dg1 = mPf.pf_fold(pString, struct1, mTemperature);

        mDDG4Ret.mD1Array[pRetIdx][mDDG4Ret.mNExps] = dg1;
        mDDG4Ret.mD2Array[pRetIdx][mDDG4Ret.mNExps] = dg2;

        ddG[mDDG4Ret.mNExps] = -1.0 * (dg1 + dg2 - dg0);
        ++mDDG4Ret.mNExps;
    }

    ddGsum = -1.0 * log_of_sum_of_exps(ddG, mDDG4Ret.mNExps);
    mDDG4Ret.mDDGSum[pRetIdx] = ddGsum;
    mDDG4Ret.mP[pRetIdx] = 1 / (1 + std::exp(ddGsum / kT));   /* Probability of being bound */

    return 0;
}

double VR16DDGWorkSpace::log_of_sum_of_exps(std::vector<double> &pExps, int pNumExps)
{
    // compute log(exp(a_1)+exp(a_2)+...exp(a_n)) using:
    // max(a_1,a_2,..a_n) + log(1+exp(a_2-max(a_1,a_2,..a_n))+...exp(a_n-max(a_1,a_2,..a_n)))

    double exp_diff;
    double retVal;
    double max_index = 0;
    double remaining_sum = 1.0;
    double max_exp = pExps[0];
    int i;

    for (i = 1; i < pNumExps; ++i)
    {
        if (pExps[i] > max_exp)
        {
            max_exp = pExps[i];
            max_index = i;
        }
    }

    for (i = 0; i < pNumExps; ++i)
    {
        if (i != max_index)
        {
            exp_diff = pExps[i] - max_exp;
            if (exp_diff > MIN_EXP_DIFF)
            {
                remaining_sum += std::exp(exp_diff);
            }
        }
    }

    retVal = max_exp + std::log(remaining_sum);

    return retVal;
}

void VR16DDGWorkSpace::print_ddg4_ret_vals(int pRetIdx)
{
    std::cout << "dG0:    " << mDDG4Ret.mDG0[pRetIdx] << std::endl;

    std::cout << "dG1:    ";
    for (int i = 0; i < mDDG4Ret.mNExps; ++i)
    {
        std::cout << mDDG4Ret.mD1Array[pRetIdx][i];
        if (i != mDDG4Ret.mNExps - 1)
        {
            std::cout << ", ";
        }
        else
        {
            std::cout << std::endl;
        }
    }

    std::cout << "dG2:    ";
    for (int i = 0; i < mDDG4Ret.mNExps; ++i)
    {
        std::cout << mDDG4Ret.mD2Array[pRetIdx][i];
        if (i != mDDG4Ret.mNExps - 1)
        {
            std::cout << ", ";
        }
        else
        {
            std::cout << std::endl;
        }
    }

    std::cout << "ddGSum: " << mDDG4Ret.mDDGSum[pRetIdx] << std::endl;
    std::cout << "prob:   " << mDDG4Ret.mP[pRetIdx] << std::endl;
}

void VR16DDGWorkSpace::print_array_size(int pI)
{
    std::cout << pI << ": " << mDDG4Ret.mD1Array.size() << std::endl;
    std::cout << pI << ": "<< mDDG4Ret.mD1Array[0].size() << std::endl;
}

} // namespace vr16
