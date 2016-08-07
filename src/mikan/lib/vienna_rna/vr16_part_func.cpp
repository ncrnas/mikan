#include <mikan/lib/vienna_rna/include/vr16_part_func.hpp>
#include <iostream>
#include <algorithm>
#include <cfloat>

namespace vr16
{

//
// VR16PartFunc methods
//
void VR16PartFunc::init_pf_fold(int pLen, double pTemperature)
{
    if (pLen < 1)
    {
        std::cerr << "Error: init_pf_fold - length must be greater 0." << std::endl;
    }

    resize_arrays((unsigned) pLen);
    mPFParams.scale_pf_params((unsigned) pLen, pTemperature);
    mInitLength = pLen;
}

void VR16PartFunc::resize_arrays(unsigned int pLen)
{
    unsigned int arrSize;

    arrSize = (pLen + 1) * (pLen + 2) / 2;

    mQ.resize(arrSize);
    mQb.resize(arrSize);
    mQm.resize(arrSize);

    if (mStBack)
    {
        mQm1.resize(arrSize);
    }

    mPtype.resize((pLen + 1) * (pLen + 2) / 2);

    mQ1k.resize(pLen + 1);
    mQln.resize(pLen + 2);
    mQq.resize(pLen + 2);
    mQq1.resize(pLen + 2);
    mQqm.resize(pLen + 2);
    mQqm1.resize(pLen + 2);
    mPrmL0.resize(pLen + 2);
    mPrmL1.resize(pLen + 2);
    mPrml.resize(pLen + 2);

    mOpts.mIIndx.resize(pLen + 1);
    mJindx.resize(pLen + 1);
    for (unsigned int i = 1; i <= pLen; ++i)
    {
        mOpts.mIIndx[i] = ((pLen + 1 - i) * (pLen - i)) / 2 + pLen + 1;
        mJindx[i] = (i * (i - 1)) / 2;
    }
}

float VR16PartFunc::pf_fold(std::string &pString, std::string &pStructure, double pTemperature)
{
    int type_1;
    int type_2;
    int strLen;
    double temp;
    double el;

    double tmpQ;
    double Qmax = 0.0;

    double free_energy;
    double max_real;

    int idTmp = 0;
    int idQq = 0;
    int idQq1 = 1;
    int idQqm = 0;
    int idQqm1 = 1;
    int idPrmL0 = 0;
    int idPrmL1 = 1;

    max_real = DBL_MAX;

    strLen = (unsigned)pString.size();
    if (strLen > mInitLength)
    {
        init_pf_fold(strLen, pTemperature); /* (re)allocate space */
    }
    mPFParams.reset_scale();

    mS.resize((unsigned)strLen + 1);
    mS1.resize((unsigned)strLen + 1);
    mS[0] = strLen;
    for (int l = 1; l <= strLen; ++l)
    {
        mS[l] = mPairMat.encode_char(pString[l - 1], mOpts.mEnergySet);
        mS1[l] = mPairMat.get_alias(mS[l]); /* for mismatches of nostandard bases */
    }

    make_ptypes(pStructure);

    /*array initialization ; qb,qm,q
     qb,qm,q (i,j) are stored as ((n+1-i)*(n-i) div 2 + n+1-j */

    for (int d = 0; d <= TURN; ++d)
    {
        for (int i = 1; i <= strLen - d; ++i)
        {
            int j = i + d;
            int ij = mOpts.mIIndx[i] - j;
            mQ[ij] = 1.0 * mPFParams.mScale[d + 1];
            mQb[ij] = 0.0;
            mQm[ij] = 0.0;
        }
    }

    for (int i = 1; i <= strLen; ++i)
    {
        mQq[i] = 0.0;
        mQq1[i] = 0.0;
        mQqm[i] = 0.0;
        mQqm1[i] = 0.0;
        mPrmL0[i] = 0.0;
        mPrmL1[i] = 0.0;
        mPrml[i] = 0.0;
    }

    for (int j = TURN + 2; j <= strLen; ++j)
    {
        for (int i = j - TURN - 1; i >= 1; --i)
        {
            /* construction of partition function of segment i,j*/
            /*firstly that given i bound to j : qb(i,j) */
            int u = j - i - 1;
            int ij = mOpts.mIIndx[i] - j;
            double qbt1;

            type_1 = mPtype[ij];
            if (type_1 != 0)
            {
                /*hairpin contribution*/
                if (((type_1 == 3) || (type_1 == 4)) && mOpts.mNoClosingGU)
                {
                    qbt1 = 0;
                }
                else
                {
                    qbt1 = exp_hairpin_energy(u, type_1, mS1[i + 1], mS1[j - 1], i - 1, pString);
                    qbt1 *= mPFParams.mScale[u + 2];
                }

                /* interior loops with interior pair k,l */
                int max_k = j - TURN - 2;
                if (i + MAXLOOP + 1 < max_k)
                {
                    max_k = i + MAXLOOP + 1;
                }
                for (int k = i + 1; k <= max_k; ++k)
                {
                    int u1 = k - i - 1;
                    int min_l = ((k + TURN + 1) > (j - 1 - MAXLOOP + u1) ? (k + TURN + 1) : (j - 1 - MAXLOOP + u1));
                    for (int l = min_l; l < j; ++l)
                    {
                        type_2 = mPtype[mOpts.mIIndx[k] - l];
                        if (type_2 != 0)
                        {
                            type_2 = mPairMat.mRtype[type_2];
                            el = exp_loop_energy(u1, j - l - 1, type_1, type_2,
                                    mS1[i + 1], mS1[j - 1], mS1[k - 1], mS1[l + 1]);
                            qbt1 += mQb[mOpts.mIIndx[k] - l] * el * mPFParams.mScale[u1 + j - l + 1];
                        }
                    }
                }
                /*multiple stem loop contribution*/
                int ii = mOpts.mIIndx[i + 1]; /* ii-k=[i+1,k-1] */
                temp = 0.0;
                for (int k = i + 2; k <= j - 1; ++k)
                {
                    temp += mQm[ii - (k - 1)] * get_qqmx_val(idQqm1, k);
                }
                int tt = mPairMat.mRtype[type_1];
                qbt1 += temp * mPFParams.mExpMLContrib[tt][mS1[i + 1]][mS1[j - 1]] * mPFParams.mScale[2];

                mQb[ij] = qbt1;
            } /* end if (type_1!=0) */
            else
            {
                mQb[ij] = 0.0;
            }

            /* construction of qqm matrix containing final stem
             contributions to multiple loop partition function
             from segment i,j */
            set_qqmx_val(idQqm, i, get_qqmx_val(idQqm1, i) * mPFParams.mExpMLbase[1]);
            if (type_1)
            {
                qbt1 = mQb[ij] * mPFParams.mExpMLintern[type_1];
                if (i > 1)
                {
                    qbt1 *= mPFParams.mExpDangle5[type_1][mS1[i - 1]];
                }
                if (j < strLen)
                {
                    qbt1 *= mPFParams.mExpDangle3[type_1][mS1[j + 1]];
                }
                else if (type_1 > 2)
                {
                    qbt1 *= mPFParams.mExpTermAU;
                }
                set_qqmx_val(idQqm, i, get_qqmx_val(idQqm, i) + qbt1);
            }
            if (mQm1.size() != 0)
            {
                mQm1[mJindx[j] + i] = get_qqmx_val(idQqm, i); /* for stochastic backtracking */
            }

            /*construction of qm matrix containing multiple loop
             partition function contributions from segment i,j */
            temp = 0.0;
            int ii = mOpts.mIIndx[i]; /* ii-k=[i,k-1] */
            for (int k = i + 1; k <= j; ++k)
            {
                temp += (mQm[ii - (k - 1)] + mPFParams.mExpMLbase[k - i]) * get_qqmx_val(idQqm, k);
            }
            mQm[ij] = (temp + get_qqmx_val(idQqm, i));

            /*auxiliary matrix qq for cubic order q calculation below */
            qbt1 = mQb[ij];
            if (type_1)
            {
                if (i > 1)
                {
                    qbt1 *= mPFParams.mExpDangle5[type_1][mS1[i - 1]];
                }
                if (j < strLen)
                {
                    qbt1 *= mPFParams.mExpDangle3[type_1][mS1[j + 1]];
                }
                else if (type_1 > 2)
                {
                    qbt1 *= mPFParams.mExpTermAU;
                }
            }
            set_qqx_val(idQq, i, get_qqx_val(idQq1, i) * mPFParams.mScale[1] + qbt1);

            /*construction of partition function for segment i,j */
            temp = 1.0 * mPFParams.mScale[1 + j - i] + get_qqx_val(idQq, i);
            for (int k = i; k <= j - 1; ++k)
            {
                temp += mQ[ii - k] * get_qqx_val(idQq, k + 1);
            }
            mQ[ij] = temp;

            if (temp > Qmax)
            {
                Qmax = temp;
                if (Qmax > max_real / 10.0)
                {
                    std::cerr << "Warning: Q close to overflow: " << i << " " << j << " " <<  temp << std::endl;
                }
            }
            if (temp >= max_real)
            {
                std::cerr << "Error: Overflow in pf_fold while calculating q[" << i << "," << j << "]" << std::endl;
                std::cerr << "use larger pf_scale" << std::endl;
            }
        }
        idTmp = idQq1;
        idQq1 = idQq;
        idQq = idTmp;

        idTmp = idQqm1;
        idQqm1 = idQqm;
        idQqm = idTmp;

    }

    if (mOpts.mBacktrackType == 'C')
    {
        tmpQ = mQb[mOpts.mIIndx[1] - strLen];
    }
    else if (mOpts.mBacktrackType == 'M')
    {
        tmpQ = mQm[mOpts.mIIndx[1] - strLen];
    }
    else
    {
        tmpQ = mQ[mOpts.mIIndx[1] - strLen];
    }

    /* ensemble free energy in Kcal/mol */
    if (tmpQ <= FLT_MIN)
    {
        std::cerr << "pf_scale too large." << std::endl;
    }

    free_energy = (-1.0 * std::log(tmpQ) - strLen * std::log(mOpts.mPfScale)) *
            (pTemperature + mEn.K0) * mEn.GASCONST / 1000.0;

    /* in case we abort because of floating point errors */
    if (strLen > 1600)
    {
        std::cerr << "free energy = " << free_energy << std::endl;
    }

    /* backtracking to construct binding probabilities of pairs*/
    if (mOpts.mDoBacktrack)
    {
        backtrack(strLen, idPrmL0, idPrmL1);
        if (pStructure.size() != 0)
        {
            sprintf_bppm(strLen, pStructure);
        }
    } /* end if (do_backtrack)*/

    return (float)free_energy;
}

void VR16PartFunc::backtrack(int pStrLen, int pIdPrmL0, int pIdPrmL1)
{
    int ov = 0;
    int idTmp = 0;
    int Qmax = 0;
    double el;
    double max_real = DBL_MAX;

    for (int k = 1; k <= pStrLen; ++k)
    {
        mQ1k[k] = mQ[mOpts.mIIndx[1] - k];
        mQln[k] = mQ[mOpts.mIIndx[k] - pStrLen];
    }
    mQ1k[0] = 1.0;
    mQln[pStrLen + 1] = 1.0;

    mOpts.mProb = mQ; /* recycling */

    /* 1. exterior pair i,j and initialization of pr array */
    for (int i = 1; i <= pStrLen; ++i)
    {
        int max_j = i + TURN;
        if (pStrLen < max_j)
        {
            max_j = pStrLen;
        }
        for (int j = i; j <= max_j; ++j)
        {
            mOpts.mProb[mOpts.mIIndx[i] - j] = 0;
        }

        for (int j = i + TURN + 1; j <= pStrLen; ++j)
        {
            int ij = mOpts.mIIndx[i] - j;
            int type_1 = mPtype[ij];
            if (type_1 && (mQb[ij] > 0.0))
            {
                mOpts.mProb[ij] = mQ1k[i - 1] * mQln[j + 1] / mQ1k[pStrLen];
                if (i > 1)
                {
                    mOpts.mProb[ij] *= mPFParams.mExpDangle5[type_1][mS1[i - 1]];
                }
                if (j < pStrLen)
                {
                    mOpts.mProb[ij] *= mPFParams.mExpDangle3[type_1][mS1[j + 1]];
                }
                else if (type_1 > 2)
                {
                    mOpts.mProb[ij] *= mPFParams.mExpTermAU;
                }
            }
            else
            {
                mOpts.mProb[ij] = 0;
            }
        }
    }

    for (int l = pStrLen; l > TURN + 1; --l)
    {

        /* 2. bonding k,l as substem of 2:loop enclosed by i,j */
        for (int k = 1; k < l - TURN; k++)
        {
            int kl = mOpts.mIIndx[k] - l;
            int type_2 = mPtype[kl];
            type_2 = mPairMat.mRtype[type_2];
            if (mQb[kl] == 0)
            {
                continue;
            }

            int min_i  = k - MAXLOOP - 1;
            if (1 > min_i)
            {
                min_i = 1;
            }
            for (int i = min_i; i <= k - 1; ++i)
            {
                int max_j = l + MAXLOOP - k + i + 2;
                if (pStrLen < max_j)
                {
                    max_j = pStrLen;
                }
                for (int j = l + 1; j <= max_j; ++j)
                {
                    int ij = mOpts.mIIndx[i] - j;
                    int type_1 = mPtype[ij];
                    if ((mOpts.mProb[ij] > 0))
                    {
                        el = exp_loop_energy(k - i - 1, j - l - 1, type_1, type_2,
                                mS1[i + 1], mS1[j - 1], mS1[k - 1], mS1[l + 1]);
                        mOpts.mProb[kl] += mOpts.mProb[ij] * el * mPFParams.mScale[k - i + j - l];
                    }
                }
            }
        }
        /* 3. bonding k,l as substem of multi-loop enclosed by i,j */
        double prmMLb = 0.0;
        if (l < pStrLen)
        {
            for (int k = 2; k < l - TURN; ++k)
            {
                int i = k - 1;
                double prmt = 0.0;
                double prmt1;

                int ii = mOpts.mIIndx[i]; /* ii-j=[i,j]     */
                int ll = mOpts.mIIndx[l + 1]; /* ll-j=[l+1,j-1] */
                int tt = mPtype[ii - (l + 1)];
                tt = mPairMat.mRtype[tt];
                prmt1 = mOpts.mProb[ii - (l + 1)] * mPFParams.mExpMLclosing * mPFParams.mExpMLintern[tt]
                        * mPFParams.mExpDangle3[tt][mS1[i + 1]] * mPFParams.mExpDangle5[tt][mS1[l]];
                for (int j = l + 2; j <= pStrLen; j++)
                {
                    tt = mPtype[ii - j];
                    tt = mPairMat.mRtype[tt];
                    prmt += mOpts.mProb[ii - j] * mPFParams.mExpDangle3[tt][mS1[i + 1]]
                            * mPFParams.mExpDangle5[tt][mS1[j - 1]] * mQm[ll - (j - 1)];
                }
                int kl = mOpts.mIIndx[k] - l;
                tt = mPtype[kl];
                mPrml[i] = prmt * mPFParams.mExpMLclosing * mPFParams.mExpMLintern[tt];
                set_prmlx_val(pIdPrmL0, i, get_prmlx_val(pIdPrmL1, i) * mPFParams.mExpMLbase[1] + prmt1);

                prmMLb = prmMLb * mPFParams.mExpMLbase[1] + mPrml[i];
                /* same as:    prm_MLb = 0;
                 for (i=1; i<=k-1; i++) prm_MLb += prml[i]*expMLbase[k-i-1]; */

                mPrml[i] = mPrml[i] + get_prmlx_val(pIdPrmL0, i);

                if (mQb[kl] == 0.)
                {
                    continue;
                }

                double temp = prmMLb;

                for (i = 1; i <= k - 2; ++i)
                {
                    temp += mPrml[i] * mQm[mOpts.mIIndx[i + 1] - (k - 1)];
                }

                temp *= mPFParams.mExpMLintern[tt] * mPFParams.mScale[2];
                if (k > 1)
                {
                    temp *= mPFParams.mExpDangle5[tt][mS1[k - 1]];
                }
                if (l < pStrLen)
                {
                    temp *= mPFParams.mExpDangle3[tt][mS1[l + 1]];
                }
                mOpts.mProb[kl] += temp;

                if (mOpts.mProb[kl] > Qmax)
                {
                    Qmax = (int)mOpts.mProb[kl];
                    if (Qmax > max_real / 10.0)
                    {
                        std::cerr << "Warning: P close to overflow: " << i << " " << mOpts.mProb[kl];
                        std::cerr << " " <<  mQb[kl] << std::endl;
                    }
                }
                if (mOpts.mProb[kl] >= max_real)
                {
                    ov++;
                    mOpts.mProb[kl] = FLT_MAX;
                }

            } /* end for (k=..) */
        }

        idTmp = pIdPrmL1;
        pIdPrmL1 = pIdPrmL0;
        pIdPrmL0 = idTmp;

    } /* end for (l=..)   */

    for (int i = 1; i <= pStrLen; ++i)
    {
        for (int j = i + TURN + 1; j <= pStrLen; ++j)
        {
            int ij = mOpts.mIIndx[i] - j;
            mOpts.mProb[ij] *= mQb[ij];
        }
    }

    if (ov > 0)
    {
        std::cerr << "Warning: " << ov << " overflows occurred while backtracking;" << std::endl;
        std::cerr << "You might try a smaller pf_scale than " << mOpts.mPfScale << std::endl;
    }
}

double VR16PartFunc::get_qqx_val(int pQQId, int pIdx)
{
    double retVal = 0.0;

    if (pQQId == 0)
    {
        retVal = mQq[pIdx];
    }
    else if (pQQId == 1)
    {
        retVal = mQq1[pIdx];
    }

    return retVal;
}

void VR16PartFunc::set_qqx_val(int pQQId, int pIdx, double pNewVal)
{
    if (pQQId == 0)
    {
        mQq[pIdx] = pNewVal;
    }
    else if (pQQId == 1)
    {
        mQq1[pIdx] = pNewVal;
    }
}

double VR16PartFunc::get_qqmx_val(int pQQMId, int pIdx)
{
    double retVal = 0.0;

    if (pQQMId == 0)
    {
        retVal = mQqm[pIdx];
    }
    else if (pQQMId == 1)
    {
        retVal = mQqm1[pIdx];
    }

    return retVal;
}

void VR16PartFunc::set_qqmx_val(int pQQMId, int pIdx, double pNewVal)
{
    if (pQQMId == 0)
    {
        mQqm[pIdx] = pNewVal;
    }
    else if (pQQMId == 1)
    {
        mQqm1[pIdx] = pNewVal;
    }
}

double VR16PartFunc::get_prmlx_val(int pPrmLId, int pIdx)
{
    double retVal = 0.0;

    if (pPrmLId == 0)
    {
        retVal = mPrmL0[pIdx];
    }
    else if (pPrmLId == 1)
    {
        retVal = mPrmL1[pIdx];
    }

    return retVal;
}

void VR16PartFunc::set_prmlx_val(int pPrmLId, int pIdx, double pNewVal)
{
    if (pPrmLId == 0)
    {
        mPrmL0[pIdx] = pNewVal;
    }
    else if (pPrmLId == 1)
    {
        mPrmL1[pIdx] = pNewVal;
    }
}

void VR16PartFunc::make_ptypes(std::string &pStructure)
{
    int strLen = mS[0];

    for (int k = 1; k < strLen - TURN; k++)
    {
        for (int l = 1; l <= 2; l++)
        {
            int type_1 = 0;
            int ntype = 0;
            int otype = 0;
            int i = k;
            int j = i + TURN + l;
            if (j > strLen)
            {
                continue;
            }
            type_1 = mPairMat.mPair[mS[i]][mS[j]];
            while ((i >= 1) && (j <= strLen))
            {
                if ((i > 1) && (j < strLen))
                {
                    ntype = mPairMat.mPair[mS[i - 1]][mS[j + 1]];
                }
                if (mOpts.mNoLonelyPairs && (!otype) && (!ntype))
                {
                    type_1 = 0; /* i.j can only form isolated pairs */
                }
                mQb[mOpts.mIIndx[i] - j] = 0.;
                mPtype[mOpts.mIIndx[i] - j] = type_1;
                otype = type_1;
                type_1 = ntype;
                --i;
                ++j;
            }
        }
    }

    if (mOpts.mFoldConstrained && (pStructure.size() != 0))
    {
        int hx;
        std::vector<int> tmp_stack;
        char type_c;
        tmp_stack.resize((unsigned)strLen + 1);
        int i;
        int j;

        for (hx = 0, j = 1; j <= strLen; j++)
        {
            switch (pStructure[j - 1])
            {
            case 'x': /* can't pair */
                for (int l = 1; l < j - TURN; l++)
                {
                    mPtype[mOpts.mIIndx[l] - j] = 0;
                }
                for (int l = j + TURN + 1; l <= strLen; l++)
                {
                    mPtype[mOpts.mIIndx[j] - l] = 0;
                }
                break;
            case '(':
                tmp_stack[hx++] = j;
                for (int l = 1; l < j - TURN; l++)
                {
                    mPtype[mOpts.mIIndx[l] - j] = 0;
                }
                break;
            case '<': /* pairs upstream */
                for (int l = 1; l < j - TURN; l++)
                {
                    mPtype[mOpts.mIIndx[l] - j] = 0;
                }
                break;
            case ')':
                if (hx <= 0)
                {
                    std::cerr << "Error: Unbalanced brackets in constraints. " << pStructure << std::endl;
                }
                i = tmp_stack[--hx];
                type_c = (char)(mPtype[mOpts.mIIndx[i] - j]);
                /* don't allow pairs i<k<j<l */
                for (int k = i; k <= j; k++)
                {
                    for (int l = j; l <= strLen; l++)
                    {
                        mPtype[mOpts.mIIndx[k] - l] = 0;
                    }
                }
                /* don't allow pairs k<i<l<j */
                for (int k = 1; k <= i; k++)
                {
                    for (int l = i; l <= j; l++)
                    {
                        mPtype[mOpts.mIIndx[k] - l] = 0;
                    }
                }
                mPtype[mOpts.mIIndx[i] - j] = (type_c == 0) ? 7 : type_c;
                for (int l = j + TURN + 1; l <= strLen; l++)
                {
                    mPtype[mOpts.mIIndx[j] - l] = 0;
                }
                break;
            case '>': /* pairs downstream */
                for (int l = j + TURN + 1; l <= strLen; l++)
                {
                    mPtype[mOpts.mIIndx[j] - l] = 0;
                }
                break;
             default:
                 break;
            }
        }

        if (hx != 0)
        {
            std::cerr << "Error: Unbalanced brackets in constraints. " << pStructure << std::endl;
        }

    }
}

double VR16PartFunc::exp_hairpin_energy(int pU, int pType, int pSi1, int pSj1, int pIdx, std::string &pString)
{
    double q;
    std::string subStr6;
    std::string subStr5;

    q = mPFParams.mExpHairpin[pU];
    if ((mOpts.mTetraLoop) && (pU == 4))
    {
        subStr6 = pString.substr((unsigned)pIdx, 6);
        std::map<std::string, int>::const_iterator found = mEn.mTetraloops.find(subStr6);
        if (found != mEn.mTetraloops.end())
        {
            q *= mPFParams.mExpTetra[(*found).second];
        }
    }

    if (pU == 3)
    {
        if (mOpts.mTriLoop)
        {
            subStr5 = pString.substr((unsigned)pIdx, 5);
            std::map<std::string, int>::const_iterator found = mEn.mTriloops.find(subStr5);
            if (found != mEn.mTriloops.end())
            {
                q *= mPFParams.mExpTriloop[(*found).second];
            }
        }

        if (pType > 2)
        {
            q *= mPFParams.mExpTermAU;
        }
    }
    else
    {
        /* no mismatches for tri-loops */
        q *= mPFParams.mExpMismatchH[pType][pSi1][pSj1];
    }

    return q;
}

double VR16PartFunc::exp_loop_energy(int pU1, int pU2, int pType, int pType2, int pSi1, int pSj1, int pSp1, int pSq1)
{
    double z = 0;
    int no_close = 0;
    int absU12;

    if ((mOpts.mNoClosingGU) && ((pType2 == 3) || (pType2 == 4) || (pType == 2) || (pType == 4)))
    {
        no_close = 1;
    }

    if ((pU1 == 0) && (pU2 == 0)) /* stack */
    {
        z = mPFParams.mExpStack[pType][pType2];
    }
    else if (no_close == 0)
    {
        if ((pU1 == 0) || (pU2 == 0))
        { /* bulge */
            int u;
            u = (pU1 == 0) ? pU2 : pU1;
            z = mPFParams.mExpBulge[u];
            if (pU2 + pU1 == 1)
            {
                z *= mPFParams.mExpStack[pType][pType2];
            }
            else
            {
                if (pType > 2)
                {
                    z *= mPFParams.mExpTermAU;
                }
                if (pType2 > 2)
                {
                    z *= mPFParams.mExpTermAU;
                }
            }
        }
        else
        { /* interior loop */
            int u1u2 = pU1 + pU2;
            if (u1u2 == 2) /* size 2 is special */
            {
                z = mPFParams.mExpInt11[pType][pType2][pSi1][pSj1];
            }
            else if ((pU1 == 1) && (pU2 == 2))
            {
                z = mPFParams.mExpInt21[pType][pType2][pSi1][pSq1][pSj1];
            }
            else if ((pU1 == 2) && (pU2 == 1))
            {
                z = mPFParams.mExpInt21[pType2][pType][pSq1][pSi1][pSp1];
            }
            else if ((pU1 == 2) && (pU2 == 2))
            {
                z = mPFParams.mExpInt22[pType][pType2][pSi1][pSp1][pSq1][pSj1];
            }
            else
            {
                absU12 = pU1 - pU2;
                if (absU12 < 0)
                {
                    absU12 *= -1;
                }

                z = mPFParams.mExpInternal[u1u2] * mPFParams.mExpMismatchI[pType][pSi1][pSj1]
                      * mPFParams.mExpMismatchI[pType2][pSq1][pSp1] * mPFParams.mExpNinio[2][absU12];
            }
        }
    }

    return z;
}

void VR16PartFunc::sprintf_bppm(int pLen, std::string &pStructure)
{
    std::vector<float> P(PLEN); /* P[][0] unpaired, P[][1] upstream p, P[][2] downstream p */

    for (int j = 1; j <= pLen; ++j)
    {
        P[0] = 1.0;
        P[1] = 0.0;
        P[2] = 0.0;
        for (int i = 1; i < j; ++i)
        {
            P[2] += mOpts.mProb[mOpts.mIIndx[i] - j]; /* j is paired downstream */
            P[0] -= mOpts.mProb[mOpts.mIIndx[i] - j]; /* j is unpaired */
        }
        for (int i = j + 1; i <= pLen; ++i)
        {
            P[1] += mOpts.mProb[mOpts.mIIndx[j] - i]; /* j is paired upstream */
            P[0] -= mOpts.mProb[mOpts.mIIndx[j] - i]; /* j is unpaired */
        }
        pStructure[j - 1] = bppm_symbol(P);
    }
    pStructure[pLen] = '\0';
}

char VR16PartFunc::bppm_symbol(std::vector<float>& pP)
{
    if (pP[0] > 0.667)
    {
        return '.';
    }

    if (pP[1] > 0.667)
    {
        return '(';
    }

    if (pP[2] > 0.667)
    {
        return ')';
    }

    if ((pP[1] + pP[2]) > pP[0])
    {
        if ((pP[1] / (pP[1] + pP[2])) > 0.667)
        {
            return '{';
        }

        if ((pP[2] / (pP[1] + pP[2])) > 0.667)
        {
            return '}';
        }
        else
        {
            return '|';
        }
    }

    if (pP[0] > (pP[1] + pP[2]))
    {
        return ',';
    }

    return ':';
}

} // namespace vr16
