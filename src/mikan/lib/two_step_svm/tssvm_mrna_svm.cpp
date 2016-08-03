#include <mikan/lib/two_step_svm/include/tssvm_inst_template.hpp>  // TRNATYPE
#include <mikan/lib/two_step_svm/include/tssvm_mrna_feature.hpp>   // TSSVMRNARawFeatures
#include <mikan/lib/two_step_svm/include/tssvm_mrna_svm.hpp>       // TSSVMRNAModel, TSSVMRNAInputVector

using namespace seqan;

namespace tssvm{
//
// TSSVMRNAModel methods
//

template <class TRNAString>
void TSSVMRNAModel<TRNAString>::init_hyperplane()
{
    mHyperPlane[0] = -1.30572107595f;
    mHyperPlane[1] = 0.593744401172f;
    mHyperPlane[2] = 4.16198979669f;
    mHyperPlane[3] = 16.3787970139f;
    mHyperPlane[4] = 9.59984620853f;
    mHyperPlane[5] = 3.81289601546f;
    mHyperPlane[6] = -5.19816318534f;
    mHyperPlane[7] = -1.31499609314f;
    mHyperPlane[8] = -2.1582603992f;
    mHyperPlane[9] = -8.02015725935f;
    mHyperPlane[10] = -4.85249370022f;
    mHyperPlane[11] = -7.65730980842f;
    mHyperPlane[12] = -0.00668173334258f;
    mHyperPlane[13] = -2.95566934722f;
    mHyperPlane[14] = -12.5387352781f;
    mHyperPlane[15] = -12.0502325186f;
    mHyperPlane[16] = -4.7790370469f;
    mHyperPlane[17] = 2.79309363131f;
    mHyperPlane[18] = 4.2980051917f;
    mHyperPlane[19] = 4.46664757476f;
    mHyperPlane[20] = 5.2286212368f;
    mHyperPlane[21] = 2.34606717403f;
    mHyperPlane[22] = 1.30040788589f;
    mHyperPlane[23] = 0.889013323422f;
    mHyperPlane[24] = 2.9761142572f;
    mHyperPlane[25] = -0.188642775904f;
    mHyperPlane[26] = 0.989538955302f;
    mHyperPlane[27] = 7.82498328245f;
    mHyperPlane[28] = 0.0f;
    mHyperPlane[29] = 3.06558310165f;
    mHyperPlane[30] = 2.32548284281f;
    mHyperPlane[31] = 1.23598264889f;
    mHyperPlane[32] = 1.56331024459f;
    mHyperPlane[33] = -0.0724823473571f;
    mHyperPlane[34] = 20.5493145756f;
}

template <class TRNAString>
float TSSVMRNAModel<TRNAString>::calc_score(Eigen::VectorXf& pInput)
{
    float score;

    score = mHyperPlane.dot(pInput) + mB;

    return score;

}

//
// TSSVMRNAInputVector methods
//
template <class TRNAString>
void TSSVMRNAInputVector<TRNAString>::clear_scores()
{
    clear(mScores);
}

template <class TRNAString>
int TSSVMRNAInputVector<TRNAString>::classify(TSSVMRNARawFeatures<TRNAString> &pRNAFeatures)
{
    calc_score(pRNAFeatures);

    return 0;
}

template <class TRNAString>
int TSSVMRNAInputVector<TRNAString>::calc_score(TSSVMRNARawFeatures<TRNAString> &pRNAFeatures)
{
    TFeatSet& urlLen = pRNAFeatures.get_all_utr_len();
    TFeatSet& siteNum = pRNAFeatures.get_all_site_num();
    TFeatSet& utrLen = pRNAFeatures.get_all_tot_disc_utr_len();
    TFeatSet& seedTypeNum = pRNAFeatures.get_all_seed_type_num();
    TFeatSet& discBin =  pRNAFeatures.get_all_disc_num();
    TFeatSet& optDist = pRNAFeatures.get_all_opt_dist();
    TFeatSet& siteNumFlg = pRNAFeatures.get_all_site_num_flag();
    TFeatSet& totDisc = pRNAFeatures.get_all_to_disc();

    int k;
    resize(mScores, length(pRNAFeatures.mEffectiveRNAs));

    for (unsigned i = 0; i < length(mScores); ++i)
    {

        if (!pRNAFeatures.mEffectiveRNAs[i])
        {
            continue;
        }

        k = 0;
        mInputVec[k] = urlLen[i][0];
        ++k;
        mInputVec[k] = siteNum[i][0];
        ++k;
        mInputVec[k] = utrLen[i][0];
        ++k;
        for (unsigned j = 0; j < length(seedTypeNum[i]); ++j)
        {
            mInputVec[k] = seedTypeNum[i][j];
            ++k;
        }
        for (unsigned j = 0; j < length(discBin[i]); ++j)
        {
            mInputVec[k] = discBin[i][j];
            ++k;
        }
        for (unsigned j = 0; j < length(optDist[i]); ++j)
        {
            mInputVec[k] = optDist[i][j];
            ++k;
        }
        for (unsigned j = 0; j < length(siteNumFlg[i]); ++j)
        {
            mInputVec[k] = siteNumFlg[i][j];
            ++k;
        }
        mInputVec[k] = totDisc[i][0];

        mInputVec.normalize();
        mScores[i] = mModel.calc_score(mInputVec);

    }

    return 0;
}

template <class TRNAString>
void TSSVMRNAInputVector<TRNAString>::print_input_vector()
{
    std::cout << mInputVec.transpose();
    std::cout << std::endl;
}

// Explicit template instantiation
template class TSSVMRNAModel<TRNATYPE>;
template class TSSVMRNAInputVector<TRNATYPE>;

} // namespace tssvm
