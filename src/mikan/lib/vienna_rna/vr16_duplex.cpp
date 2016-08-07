#include <mikan/lib/vienna_rna/include/vr16_duplex.hpp>
#include <iostream>
#include <cstring>

namespace vr16
{

//
// VR16Duplex methods
//
int VR16Duplex::duplexfold(std::string &pS1, std::string &pS2, std::vector<int> &pArrayI, std::vector<int> &pArrayJ)
{
    int Emin = mParams.INF;
    int i_min = 0;
    int j_min = 0;
    int retVal;
    int dangle_energy = 0;
    int pair_type;
    int pair_type2;
    int E;
    int array_size_i = (int)pArrayI.size();
    int array_size_j = (int)pArrayJ.size();

    mBonusGiven = 0;
    mStructure.clear();

    if (mAutoFivePrimeLength == 1)
    {

        mFivePrimeLength = array_size_i;
        if (mDebugMode)
        {
            std::cerr << "Five prime length set to array size = " << mFivePrimeLength << std::endl;
        }
    }

    if (mDebugMode)
    {
        std::cerr << "array_size_i = "<< array_size_i << ", array_size_j = " << array_size_j << std::endl;

        std::cerr << "array_i = [ ";
        for (int i = 0; i < array_size_i; ++i)
        {
            std::cerr << pArrayI[i] << " ";
        }
        std::cerr << "]" << std::endl;
        std::cerr << "array_j = [ ";
        for (int i = 0; i < array_size_j; ++i)
        {
            std::cerr << pArrayJ[i] << " ";
        }
        std::cerr << "]" << std::endl;
    }

    mN1 = (int)pS1.size();
    mN2 = (int)pS2.size();

    mC.resize((unsigned)mN1 + 1);

    for (int i = 0; i <= mN1; ++i)
    {
        mC[i].resize((unsigned)mN2 + 1);
    }

    encode_seq(pS1, pS2);

    for (int i = 1; i <= mN1; ++i)
    {

        for (int j = mN2; j > 0; --j)
        {
            pair_type = mPairMat.mPair[mS1[i]][mS2[j]];
            mC[i][j] = pair_type ? mParams.mDuplexInit : mParams.INF;

            if (pair_type == 0)
            {
                if ((i <= array_size_i) && (pArrayI[i - 1] == mN2 - j + 1))
                {
                    std::cerr << "ERROR: Wrong constraint given at " << i << "," << j << std::endl;
                }

                continue;
            }

            if (array_size_i)
            {
                if ((i <= array_size_i) && (pArrayI[i - 1] == mN2 - j + 1))
                {
                    mC[i][j] -= BONUS_SIZE;
                    mBonusGiven = BONUS_SIZE;

                    if (mDebugMode)
                    {
                        std::cerr << "Gave bonus to " << i << "," << j << std::endl;
                    }
                }

                if ((i <= array_size_i && (pArrayI[i - 1] != mN2 - j + 1))
                        || (mN2 - j + 1 <= array_size_j && (pArrayJ[mN2 - j] != i)))
                {
                    if (mDebugMode > 1)
                    {
                        std::cerr << "> Skipping: " << i << "," << j << std::endl;
                    }

                    mC[i][j] = mParams.INF;
                    continue;
                }
            }

            if (i > 1)
            {
                mC[i][j] += mParams.mDangle5[pair_type][mSS1[i - 1]];
            }
            if (j < mN2)
            {
                mC[i][j] += mParams.mDangle3[pair_type][mSS2[j + 1]];
            }
            if (pair_type > 2)
            {
                mC[i][j] += mParams.mTerminalAU;
            }

            for (int k = i - 1; k > 0 && k > i - mParams.MAXLOOP - 2; --k)
            {
                for (int l = j + 1; l <= mN2; ++l)
                {
                    if (i - k + l - j - 2 >  mParams.MAXLOOP)
                    {
                        break;
                    }

                    pair_type2 = mPairMat.mPair[mS1[k]][mS2[l]];
                    if (!pair_type2)
                    {
                        continue;
                    }

                    if (mC[i][j] < mC[k][l] +  mPpIL.mMinIL[i - k - 1][l - j - 1])
                    {
                        continue;
                    }

                    E = mPpIL.loop_energy(i - k - 1, l - j - 1, pair_type2, mPairMat.mRtype[pair_type],
                            mSS1[k + 1], mSS2[l - 1], mSS1[i - 1], mSS2[j + 1]);

                    mC[i][j] = ((mC[i][j] < (mC[k][l] + E)) ? mC[i][j] : (mC[k][l] + E));
                }
            }

            E = mC[i][j];

            if (i < mN1)
            {
                E += mParams.mDangle3[mPairMat.mRtype[pair_type]][mSS1[i + 1]];
            }
            if (j > 1)
            {
                E += mParams.mDangle5[mPairMat.mRtype[pair_type]][mSS2[j - 1]];
            }
            if (pair_type > 2)
            {
                E += mParams.mTerminalAU;
            }

            if (E < Emin)
            {
                Emin = E;
                i_min = i;
                j_min = j;
            }
        }
    }

    if (mDebugMode)
    {
        /* Print dynamic programming matrix */
        for (int j = 1; j <= mN2; ++j)
        {
            std::cerr << "\t" << j << ",";
        }
        std::cerr << std::endl;

        for (int i = 1; i <= mN1; i++)
        {
            std::cerr << i;

            for (int j = 1; j <= mN2; j++)
            {
                if (mC[i][j] == mParams.INF)
                {
                    std::cerr << "\tINF";
                }
                else
                {
                    std::cerr << "\t" << mC[i][j];
                }
            }
            std::cerr << std::endl;
        }
    }

    if (mOpts.mDoBacktrack)
    {
        retVal = backtrack(i_min, j_min, pArrayI, pArrayJ);
        if (retVal != 0)
        {
            return 1;
        }
    }

    if (i_min < mN1)
    {
        ++i_min;
    }
    if (j_min > 1)
    {
        --j_min;
    }

    mI = i_min;
    mJ = j_min;

    if (array_size_i)
    {
        Emin += mBonusGiven;
    }

    mEnergy = (float) (Emin) / 100.0f;

    if (mFivePrimeLength)
    {
        dangle_energy = Emin - (mE3p + mE5p);

        if (mDebugMode)
        {
            std::cerr << "Dangle energy: " << dangle_energy << std::endl;
        }

        if (i_min > mFivePrimeLength)
        {
            mE3p += dangle_energy;
        }
        else
        {
            mE5p += dangle_energy;
        }
    }

    if (!mDelayFree)
    {
        for (int i = 0; i <= mN1; ++i)
        {
            mC[i].clear();
        }
        mC.clear();
        mS1.clear();
        mS2.clear();
        mSS1.clear();
        mSS2.clear();
    }

    return 0;
}

void VR16Duplex::encode_seq(std::string &pS1, std::string& pS2)
{
    unsigned l;

    l = (unsigned)pS1.size();
    mS1.resize(l + 1, 0);
    mSS1.resize(l + 1, 0);

    /* SS1 exists only for the special X K and I bases and energy_set!=0 */
    for (unsigned i = 1; i <= l; ++i)
    { /* make numerical encoding of sequence */
        mS1[i] = mPairMat.encode_char(pS1[i - 1], mOpts.mEnergySet);
        mSS1[i] = mPairMat.get_alias(mS1[i]); /* for mismatches of nostandard bases */
    }

    l = (unsigned)pS2.size();
    mS2.resize(l + 1, 0);
    mSS2.resize(l + 1, 0);

    /* SS2 exists only for the special X K and I bases and energy_set!=0 */
    for (unsigned i = 1; i <= l; ++i)
    { /* make numerical encoding of sequence */
        mS2[i] = mPairMat.encode_char(pS2[i - 1], mOpts.mEnergySet);
        mSS2[i] = mPairMat.get_alias(mS2[i]); /* for mismatches of nostandard bases */
    }
}

int VR16Duplex::backtrack(int pI, int pJ, std::vector<int> &pArrayI, std::vector<int> &pArrayJ)
{

    /* backtrack structure going backwards from i, and forwards from j
     return structure in bracket notation with & as separator */

    int pair_type;
    int pair_type2;
    int E;
    int LE;
    int traced;
    int i0;
    int j0;
    std::string st1(mN1 + 1, 0);
    std::string st2(mN2 + 1, 0);
    int array_size_i = (int)pArrayI.size();
    int array_size_j = (int)pArrayJ.size();

    int first_time_five_prime = 1;
    int e_min;

    e_min = mC[pI][pJ];

    if (mDebugMode)
    {
        std::cerr << "Starting backtrack at " << pI << ","<< pJ << " with " << mC[pI][pJ] << std::endl;
    }

    i0 = (((pI + 1) < mN1) ? (pI + 1) : mN1);
    j0 = (((pJ - 1) > 1) ? (pJ - 1) : 1);

    while (pI > 0 && pJ <= mN2)
    {
        E = mC[pI][pJ];
        traced = 0;
        st1[pI - 1] = '(';
        st2[pJ - 1] = ')';

        pair_type = mPairMat.mPair[mS1[pI]][mS2[pJ]];

        if (pair_type == 0)
        {
            std::cerr << "ERROR: Backtrack failed in fold duplex." << std::endl;
            return 1;
        }

        for (int k = pI - 1; k > 0 && k > pI - mParams.MAXLOOP - 2; --k)
        {
            for (int l = pJ + 1; l <= mN2; ++l)
            {
                if (pI - k + l - pJ - 2 >  mParams.MAXLOOP)
                {
                    break;
                }

                pair_type2 = mPairMat.mPair[mS1[k]][mS2[l]];
                if (!pair_type2)
                {
                    continue;
                }

                if (array_size_i)
                {
                    if ((k <= array_size_i && (pArrayI[k - 1] != mN2 - l + 1))
                            || (mN2 - l + 1 <= array_size_j && (pArrayJ[mN2 - l] != k)))
                    {
                        continue;
                    }
                }

//                if (E < mC[k][l] +  mPpIL.mMinIL[pI - k - 1][l - pJ - 1])
//                {
//                    continue;
//                }

                LE = mPpIL.loop_energy(pI - k - 1, l - pJ - 1, pair_type2, mPairMat.mRtype[pair_type],
                        mSS1[k + 1], mSS2[l - 1], mSS1[pI - 1], mSS2[pJ + 1]);

                if (E == mC[k][l] + LE)
                {
                    traced = 1;
                    pI = k;
                    pJ = l;
                    if (mDebugMode)
                    {
                        std::cerr << "Traced: k=" << k << ", l=" << l << ", E=" <<  E << std::endl;
                    }

                    if (mFivePrimeLength)
                    {
                        if (k < mFivePrimeLength && first_time_five_prime)
                        {
                            first_time_five_prime = 0;
                            mE3p = e_min - E;
                            mE5p = E + (mBonusGiven * (array_size_i > 0));

                            if (mDebugMode)
                            {
                                std::cerr << "Reached five prime border " << mFivePrimeLength<< " at k=" << k;
                                std::cerr << ". Thus, E3p = " << mE3p << ", E5p = " << mE5p << std::endl;
                            }
                        }
                    }

                    break;
                }
            }

            if (traced)
            {
                break;
            }
        }

        if (!traced)
        {
            if (pI > 1)
            {
                E -= mParams.mDangle5[pair_type][mSS1[pI - 1]];
            }
            if (pJ < mN2)
            {
                E -= mParams.mDangle3[pair_type][mSS2[pJ + 1]];
            }
            if (pair_type > 2)
            {
                E -= mParams.mTerminalAU;
            }

            if (E != mParams.mDuplexInit - (mBonusGiven * (array_size_i > 0)))
            {
                std::cerr << "Error: Backtrack failed in fold duplex.";
                std::cerr << " E = " << E;
                std::cerr << ", A potential problem with impossible restrictions with -i, -j." << std::endl;
                break;
            }
            else
            {
                break;
            }
        }
    }

    if (array_size_i)
    {
        E += mBonusGiven;
    }

    if (pI > 1)
    {
        --pI;
    }
    if (pJ < mN2)
    {
        ++pJ;
    }

    for (int k = pI; k <= i0; ++k)
    {
        if (!st1[k - 1])
        {
            st1[k - 1] = '.';
        }
    }

    for (int k = j0; k <= pJ; ++k)
    {
        if (!st2[k - 1])
        {
            st2[k - 1] = '.';
        }
    }

    mStructure.resize((unsigned)(i0 - pI + 1 + pJ - j0 + 1 + 1), 0);

    int idx = 0;
    for (int i = pI - 1; st1[i] != 0; ++i)
    {
        mStructure[idx] = st1[i];
        ++idx;
    }
    mL1 = idx;
    mStructure[idx] = '&';
    ++idx;

    for (int i = j0 - 1; st2[i] != 0; ++i)
    {
        mStructure[idx] = st2[i];
        ++idx;
    }

    return 0;
}

} // namespace vr16
