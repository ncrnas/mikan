#include <vr16_params.hpp>
#include <algorithm>                  // std::min
#include <iostream>
#include <cmath>
#include <vector>

namespace vr16
{

//
// VR16Params methods
//
void VR16Params::init_parameters(double pTemperature)
{
    double tempf = ((pTemperature + K0) / mEn.mTmeasure);

    for (unsigned i = 0; i < 31; ++i)
    {
        mHairpin[i] = (int) mEn.mHairpin37[i] * (tempf);
    }

    int imax = 0;
    for (int i = 0; i <= std::min(30, (int)MAXLOOP); ++i)
    {
        mBulge[i] = (int) mEn.mBulge37[i] * tempf;
        mInternalLoop[i] = (int) mEn.mInternalLoop37[i] * tempf;
        ++imax;
    }

    mlxc = mEn.mlxc37 * tempf;
    for (int i = imax; i <= MAXLOOP; ++i)
    {
        mBulge[i] = mBulge[30] + (int) (mlxc * std::log((double) (i) / 30.0));
        mInternalLoop[i] = mInternalLoop[30] + (int) (mlxc * std::log((double) (i) / 30.0));
    }

    for (int i = 0; i < 5; ++i)
    {
        mFNinio[i] = (int) mEn.mFNinio37[i] * tempf;
    }

    for (unsigned i = 0; i < mEn.mTetraloops.size(); ++i)
    {
        mTetraEnergy[i] = mEn.mTetraEnth37
                - (mEn.mTetraEnth37 - mEn.mTetraEnergy37[i]) * tempf;
    }
    for (unsigned i = 0; i < mEn.mTriloops.size(); ++i)
    {
        mTriloopE[i] = mEn.mTriloopE37[i];
    }

    mMLBase = mEn.mMLBase37 * tempf;
    for (int i = 0; i <= NBPAIRS; ++i)
    { /* includes AU penalty */
        mMLintern[i] = mEn.mMLIntern37 * tempf;
        mMLintern[i] += (i > 2) ? mEn.mTerminalAU : 0;
    }
    mMLclosing = mEn.mMLClosing37 * tempf;

    mTerminalAU = mEn.mTerminalAU;

    mDuplexInit = mEn.mDuplexInit * tempf;

    /* stacks    G(T) = H - [H - G(T0)]*T/T0 */
    for (int i = 0; i <= NBPAIRS; ++i)
    {
        for (int j = 0; j <= NBPAIRS; ++j)
        {
            mStack[i][j] = mEn.mEnthalpies[i][j] - (mEn.mEnthalpies[i][j] - mEn.mStack37[i][j]) * tempf;
        }
    }

    /* mismatches */
    for (int i = 0; i <= NBPAIRS; ++i)
    {
        for (int j = 0; j < 5; ++j)
        {
            for (int k = 0; k < 5; ++k)
            {
                mMismatchI[i][j][k] = mEn.mMismH[i][j][k] - (mEn.mMismH[i][j][k] - mEn.mMismatchI37[i][j][k]) * tempf;
                mMismatchH[i][j][k] = mEn.mMismH[i][j][k] - (mEn.mMismH[i][j][k] - mEn.mMismatchH37[i][j][k]) * tempf;
                mMismatchM[i][j][k] = mEn.mMismH[i][j][k] - (mEn.mMismH[i][j][k] - mEn.mMismatchM37[i][j][k]) * tempf;
            }
        }
    }

    /* dangles */
    int dd;
    for (int i = 0; i <= NBPAIRS; ++i)
    {
        for (int j = 0; j < 5; ++j)
        {
            dd = mEn.mDangle5_H[i][j] - (mEn.mDangle5_H[i][j] - mEn.mDangle5_37[i][j]) * tempf;
            mDangle5[i][j] = (dd > 0) ? 0 : dd; /* must be <= 0 */

            dd = mEn.mDangle3_H[i][j] - (mEn.mDangle3_H[i][j] - mEn.mDangle3_37[i][j]) * tempf;
            mDangle3[i][j] = (dd > 0) ? 0 : dd; /* must be <= 0 */
        }
    }

    /* interior 1x1 loops */
    for (int i = 0; i <= NBPAIRS; ++i)
    {
        for (int j = 0; j <= NBPAIRS; ++j)
        {
            for (int k = 0; k < 5; ++k)
            {
                for (int l = 0; l < 5; ++l)
                {
                    mInt11[i][j][k][l] = mEn.mInt11_H[i][j][k][l] -
                            (mEn.mInt11_H[i][j][k][l] - mEn.mInt11_37[i][j][k][l]) * tempf;
                }
            }
        }
    }

    /* interior 2x1 loops */
    for (int i = 0; i <= NBPAIRS; ++i)
    {
        for (int j = 0; j <= NBPAIRS; ++j)
        {
            for (int k = 0; k < 5; ++k)
            {
                for (int l = 0; l < 5; ++l)
                {
                    for (int m = 0; m < 5; ++m)
                    {
                        mInt21[i][j][k][l][m] = mEn.mInt21_H[i][j][k][l][m] - (mEn.mInt21_H[i][j][k][l][m]
                                                 - mEn.mInt21_37[i][j][k][l][m]) * tempf;
                    }
                }
            }
        }
    }

    /* interior 2x2 loops */
    for (int i = 0; i <= NBPAIRS; ++i)
    {
        for (int j = 0; j <= NBPAIRS; ++j)
        {
            for (int k = 0; k < 5; ++k)
            {
                for (int l = 0; l < 5; ++l)
                {
                    for (int m = 0; m < 5; ++m)
                    {
                        for (int n = 0; n < 5; ++n)
                        {
                            mInt22[i][j][k][l][m][n] = mEn.mInt22_H[i][j][k][l][m][n] - (mEn.mInt22_H[i][j][k][l][m][n]
                                                        - mEn.mInt22_37[i][j][k][l][m][n]) * tempf;
                        }
                    }
                }
            }
        }
    }

    mTemperature = pTemperature;
    ++mId;

}

//
// VR16ParamIL methods
//
void VR16ParamIL::init_parameters()
{
    int lp;

    for (int i = 0; i < mParams.MAXLOOP + 1; ++i)
    {
        for (int j = 0; j < mParams.MAXLOOP + 1; ++j)
        {
            mMinIL[i][j] = INF;
            for (int k = 0; k < mParams.NBPAIRS + 1; ++k)
            {
                for (int l = 0; l < mParams.NBPAIRS + 1; ++l)
                {
                    for (int m = 0; m < 5; ++m)
                    {
                        for (int n = 0; n < 5; ++n)
                        {
                            for (int o = 0; o < 5; ++o)
                            {
                                for (int p = 0; p < 5; ++p)
                                {
                                    lp = loop_energy(i, j, k, l, m, n, o, p);
                                    if (lp < mMinIL[i][j])
                                    {
                                        mMinIL[i][j] = lp;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }

}

int VR16ParamIL::loop_energy(int pN1, int pN2, int pType1, int pType2, int pSi1, int pSj1, int pSp1, int pSq1)
{

    /* compute energy of degree 2 loop (stack bulge or interior) */
    int nl, ns, energy;

    if (pN1 > pN2)
    {
        nl = pN1;
        ns = pN2;
    }
    else
    {
        nl = pN2;
        ns = pN1;
    }

    if (nl == 0)
    {
        return mParams.mStack[pType1][pType2]; /* stack */
    }
    else if (ns == 0)
    { /* bulge */
        energy = (nl <= mParams.MAXLOOP) ? mParams.mBulge[nl] :
                (mParams.mBulge[30] + (int) (mParams.mlxc * std::log(nl / 30.0)));
        if (nl == 1)
        {
            energy += mParams.mStack[pType1][pType2];
        }
        else
        {
            if (pType1 > 2)
            {
                energy += mParams.mTerminalAU;
            }
            if (pType2 > 2)
            {
                energy += mParams.mTerminalAU;
            }
        }
        return energy;
    }
    else
    { /* interior loop */
        if (ns == 1)
        {
            if (nl == 1) /* 1x1 loop */
            {
                return mParams.mInt11[pType1][pType2][pSi1][pSj1];
            }
            if (nl == 2)
            { /* 2x1 loop */
                if (pN1 == 1)
                {
                    energy = mParams.mInt21[pType1][pType2][pSi1][pSq1][pSj1];
                }
                else
                {
                    energy = mParams.mInt21[pType2][pType1][pSq1][pSi1][pSp1];
                }
                return energy;
            }
        }
        else if (pN1 == 2 && pN2 == 2) /* 2x2 loop */
        {
            return mParams.mInt22[pType1][pType2][pSi1][pSp1][pSq1][pSj1];
        }

        int n1n2 = pN1 + pN2;
        if (n1n2 <= mParams.MAXLOOP)
        {
            energy =  mParams.mInternalLoop[n1n2];
        }
        else
        {
            energy = mParams.mInternalLoop[30] + (int) (mParams.mlxc * std::log((n1n2) / 30.0));
        }

        int energy_t = (nl - ns) * mParams.mFNinio[2];
        if (mEn.MAX_NINIO < energy_t)
        {
            energy_t = mEn.MAX_NINIO;
        }
        energy += energy_t + mParams.mMismatchI[pType1][pSi1][pSj1] + mParams.mMismatchI[pType2][pSq1][pSp1];

    }

    return energy;
}

//
// VR16PartFunc methods
//
double VR16PFParams::pf_smooth(double pX)
{
    return ((pX / PF_SCALE_10) < -1.2283697) ? 0.0 :
            (((pX / PF_SCALE_10) > 0.8660254) ? pX : 0.38490018 * PF_SCALE_10
                    * (std::sin((pX / PF_SCALE_10) - 0.34242663) + 1.0)
                    * (std::sin((pX / PF_SCALE_10) - 0.34242663) + 1.0));
}

void VR16PFParams::scale_pf_params(unsigned int pLen, double pTemperature)
{
    /* scale energy parameters and pre-calculate Boltzmann weights */
//    unsigned int i, j, k, l;
    double kT;
    double TT;
    double GT;

    if (pLen <= (unsigned) mInitLength)
    {
        return;
    }
    mInitLength = pLen;
    mInitTemp = pTemperature;

    mExpHairpin.resize(pLen + 1);
    mExpMLbase.resize(pLen + 1);
    mScale.resize(pLen + 1);

    kT = (pTemperature + mEn.K0) * mEn.GASCONST; /* kT in cal/mol  */
    TT = (pTemperature + mEn.K0) / (mEn.mTmeasure);

    /* scaling factors (to avoid overflows) */
    if (mOpts.mPfScale == -1)
    { /* mean energy for random sequences: 184.3*length cal */
        mOpts.mPfScale = std::exp(-1.0 * (-185.0 + (pTemperature - 37.0) * 7.27) / kT);
        if (mOpts.mPfScale < 1)
        {
            mOpts.mPfScale = 1;
        }
    }

    /* loop energies: hairpins, bulges, interior, mulit-loops */
    for (unsigned i = 0; i <= (unsigned)std::min(30, (int)pLen); ++i)
    {
        GT = mEn.mHairpin37[i] * TT;
        mExpHairpin[i] = std::exp(-10.0 * GT / kT);
    }

    for (unsigned i = 0; i <= (unsigned)std::min(30, (int)MAXLOOP); ++i)
    {
        GT = mEn.mBulge37[i] * TT;
        mExpBulge[i] = std::exp(-10.0 * GT / kT);
        GT = mEn.mInternalLoop37[i] * TT;
        mExpInternal[i] = std::exp(-10.0 * GT / kT);
    }

    /* special case of size 2 interior loops (single mismatch) */
    if (mOpts.mJamesRule)
    {
        mExpInternal[2] = std::exp(-80.0 * 10.0 / kT);
    }

    double lxcPf = mEn.mlxc37 * TT;
    for (unsigned i = 31; i < pLen; ++i)
    {
        GT = mEn.mHairpin37[30] * TT + (lxcPf * std::log(i / 30.0));
        mExpHairpin[i] = std::exp(-10.0 * GT / kT);
    }

    for (unsigned i = 31; i <= (unsigned)MAXLOOP; ++i)
    {
        GT = mEn.mBulge37[30] * TT + (lxcPf * std::log(i / 30.0));
        mExpBulge[i] = std::exp(-GT * 10.0 / kT);
        GT = mEn.mInternalLoop37[30] * TT + (lxcPf * std::log(i / 30.0));
        mExpInternal[i] = std::exp(-GT * 10.0 / kT);
    }

    for (unsigned i = 0; i < 5; ++i)
    {
        GT = mEn.mFNinio37[i] * TT;
        for (unsigned j = 0; j <= (unsigned)MAXLOOP; ++j)
        {
            mExpNinio[i][j] = std::exp(-10.0 * std::min(mEn.MAX_NINIO, (int)(j * GT)) / kT);
        }
    }

    for (unsigned i = 0; i < mEn.mTetraloops.size(); ++i)
    {
        GT = mEn.mTetraEnth37 - (mEn.mTetraEnth37 - mEn.mTetraEnergy37[i]) * TT;
        mExpTetra[i] = std::exp(-10.0 * GT / kT);
    }

    for (unsigned i = 0; i < mEn.mTriloops.size(); ++i)
    {
        mExpTriloop[i] = std::exp(-10.0 * mEn.mTriloopE37[i] / kT);
    }

    GT = mEn.mMLClosing37 * TT;
    mExpMLclosing = std::exp(-10.0 * GT / kT);

    for (unsigned i = 0; i <= (unsigned)NBPAIRS; ++i)
    { /* includes AU penalty */
        GT = mEn.mMLIntern37 * TT;
        /* if (i>2) GT += TerminalAU; */
        mExpMLintern[i] = std::exp(-10.0 * GT / kT);
    }
    mExpTermAU = std::exp(-10.0 * mEn.mTerminalAU / kT);

    /* if dangles==0 just set their energy to 0,
     don't let dangle energies become > 0 (at large temps),
     but make sure go smoothly to 0                        */
    for (unsigned i = 0; i <= (unsigned)NBPAIRS; ++i)
    {
        for (unsigned j = 0; j <= 4; ++j)
        {
            if (mOpts.mDangles)
            {
                /*------------------------------------------------------------------------*/
                /* dangling ends should never be destabilizing, i.e. expdangle>=1         */
                /* specific heat needs smooth function (2nd derivative)                   */
                /* we use a*(sin(x+b)+1)^2, with a=2/(3*sqrt(3)), b=Pi/6-sqrt(3)/2,       */
                /* in the interval b<x<sqrt(3)/2                                          */
                GT = mEn.mDangle5_H[i][j] - (mEn.mDangle5_H[i][j] - mEn.mDangle5_37[i][j]) * TT;
                mExpDangle5[i][j] = std::exp(pf_smooth(-1.0 * GT) * 10.0 / kT);
                GT = mEn.mDangle3_H[i][j]- (mEn.mDangle3_H[i][j] - mEn.mDangle3_37[i][j]) * TT;
                mExpDangle3[i][j] = std::exp(pf_smooth(-1.0 * GT) * 10.0 / kT);
            }
            else
            {
                mExpDangle3[i][j] = mExpDangle5[i][j] = 1;
            }

            if (i > 2) /* add TermAU penalty into dangle3 */
            {
                mExpDangle3[i][j] *= mExpTermAU;
            }
        }
    }

    /* stacking energies */
    for (unsigned i = 0; i <= (unsigned)NBPAIRS; ++i)
    {
        for (unsigned j = 0; j <= (unsigned)NBPAIRS; ++j)
        {
            GT = mEn.mEnthalpies[i][j] - (mEn.mEnthalpies[i][j] - mEn.mStack37[i][j]) * TT;
            mExpStack[i][j] = std::exp(-10.0 * GT / kT);
        }
    }

    /* mismatch energies */
    for (unsigned i = 0; i <= (unsigned)NBPAIRS; ++i)
    {
        for (unsigned j = 0; j < 5; ++j)
        {
            for (unsigned k = 0; k < 5; ++k)
            {
                GT = mEn.mMismH[i][j][k] - (mEn.mMismH[i][j][k] - mEn.mMismatchI37[i][j][k]) * TT;
                mExpMismatchI[i][j][k] = std::exp(-10.0 * GT / kT);

                GT = mEn.mMismH[i][j][k] - (mEn.mMismH[i][j][k] - mEn.mMismatchH37[i][j][k]) * TT;
                mExpMismatchH[i][j][k] = std::exp(-10.0 * GT / kT);

                GT = mEn.mMismH[i][j][k] - (mEn.mMismH[i][j][k] - mEn.mMismatchM37[i][j][k]) * TT;
                mExpMismatchM[i][j][k] = std::exp(-10.0 * GT / kT);
            }
        }
    }

    /* interior lops of length 2 */
    for (unsigned i = 0; i <= (unsigned)NBPAIRS; ++i)
    {
        for (unsigned j = 0; j <= (unsigned)NBPAIRS; ++j)
        {
            for (unsigned k = 0; k < 5; ++k)
            {
                for (unsigned l = 0; l < 5; ++l)
                {
                    GT = mEn.mInt11_H[i][j][k][l] - (mEn.mInt11_H[i][j][k][l] - mEn.mInt11_37[i][j][k][l]) * TT;
                    mExpInt11[i][j][k][l] = std::exp(-10.0 * GT / kT);
                }
            }
        }
    }

    /* interior 2x1 loops */
    for (unsigned i = 0; i <= (unsigned)NBPAIRS; ++i)
    {
        for (unsigned j = 0; j <= (unsigned)NBPAIRS; ++j)
        {
            for (unsigned k = 0; k < 5; ++k)
            {
                for (unsigned l = 0; l < 5; ++l)
                {
                    for (unsigned m = 0; m < 5; ++m)
                    {
                        GT = mEn.mInt21_H[i][j][k][l][m] - (mEn.mInt21_H[i][j][k][l][m]
                              - mEn.mInt21_37[i][j][k][l][m]) * TT;
                        mExpInt21[i][j][k][l][m] = std::exp(-10.0 * GT / kT);
                    }
                }
            }
        }
    }

    /* interior 2x2 loops */
    for (unsigned i = 0; i <= (unsigned)NBPAIRS; ++i)
    {
        for (unsigned j = 0; j <= (unsigned)NBPAIRS; ++j)
        {
            for (unsigned k = 0; k < 5; ++k)
            {
                for (unsigned l = 0; l < 5; ++l)
                {
                    for (unsigned m = 0; m < 5; ++m)
                    {
                        for (unsigned n = 0; n < 5; ++n)
                        {
                            GT = mEn.mInt22_H[i][j][k][l][m][n]
                                 - (mEn.mInt22_H[i][j][k][l][m][n]- mEn.mInt22_37[i][j][k][l][m][n]) * TT;
                            mExpInt22[i][j][k][l][m][n] = std::exp(-10.0 * GT / kT);
                        }
                    }
                }
            }
        }
    }

    for (unsigned i = 0; i <= (unsigned)NBPAIRS; ++i)
    {
        for (unsigned j = 0; j <= 4; ++j)
        {
            for (unsigned k = 0; k <= 4; ++k)
            {
                mExpMLContrib[i][j][k] = mExpMLclosing * mExpMLintern[i] * mExpDangle3[i][j] * mExpDangle5[i][k];
            }
        }
    }


}

void VR16PFParams::reset_scale()
{
    double kT;
    double TT;
    double GT;

    kT = (mInitTemp + mEn.K0) * mEn.GASCONST;         /* kT in cal/mol  */
    TT = (mInitTemp + mEn.K0) / (mEn.mTmeasure);
    GT = mEn.mMLBase37 * TT;

    mScale[0] = 1.0;
    for (unsigned i = 1; i <= mScale.size(); ++i)
    {
        mScale[i] = mScale[i - 1] / mOpts.mPfScale;
    }

    for (unsigned i = 0; i < mExpMLbase.size(); ++i)
    {
        mExpMLbase[i] = std::exp(-10.0 * i * GT / kT) * mScale[i];
    }

}

} // namespace vr16
