#include <mk_inst_template.hpp>     // TRNATYPE
#include <tssvm_site_feature.hpp>   // TSSVMRawFeatures
#include <tssvm_site_svm.hpp>       // TSSVMSiteModel, TSSVMSiteInputVector
#include <tssvm_site_svm_alpha.hpp> // init_alpha_vector
#include <tssvm_site_svm_sv.hpp>    // TSSVMSiteModelSV

using namespace seqan;
using namespace mikan;

namespace tssvm {
//
// TSSVMSiteModel methods
//

template<class TRNAString>
int TSSVMSiteModel<TRNAString>::init_model() {
    int retVal;

    retVal = init_alpha();
    if (retVal != 0) {
        return retVal;
    }

    retVal = init_sv();
    if (retVal != 0) {
        return retVal;
    }

    return 0;
}

template<class TRNAString>
int TSSVMSiteModel<TRNAString>::init_alpha() {
    return init_alpha_vector(mAlphas);
}

template<class TRNAString>
int TSSVMSiteModel<TRNAString>::init_sv() {
    TSSVMSiteModelSV sitesv(mSVs);
    sitesv.init_sv_matix();

    return 0;
}

template<class TRNAString>
float TSSVMSiteModel<TRNAString>::calc_score(Eigen::VectorXf &pInput) {
    float score;

    mMatProd = mSVs * pInput;

    for (int i = 0; i < MODEL_M; ++i) {
        mMatProd(i) = std::pow((float) mMatProd(i), K_POLY_DEGREE);
    }

    score = mAlphas.dot(mMatProd) + mB;

    return score;

}

//
// TSSVMSiteInputVector methods
//
template<class TRNAString>
void TSSVMSiteInputVector<TRNAString>::clear_scores() {
    clear(mScores);
}

template<class TRNAString>
int TSSVMSiteInputVector<TRNAString>::classify(TSSVMRawFeatures<TRNAString> &pSiteFeatures) {
    calc_score(pSiteFeatures);

    return 0;
}

template<class TRNAString>
int TSSVMSiteInputVector<TRNAString>::calc_score(TSSVMRawFeatures<TRNAString> &pSiteFeatures) {
    TFeatSet &seedTypes = pSiteFeatures.get_all_seed_type();
    TFeatSet &similarities = pSiteFeatures.get_all_similarities();
    TFeatSet &auRichUp = pSiteFeatures.get_all_au_rich_up();
    TFeatSet &auRichDown = pSiteFeatures.get_all_au_rich_down();
    TFeatSet &sitePos = pSiteFeatures.get_all_site_pos();
    TFeatSet &seqMatch = pSiteFeatures.get_all_seq_match();
    TFeatSet &a1Match = pSiteFeatures.get_all_a1_match();

    int k;
    resize(mScores, length(pSiteFeatures.mEffectiveSites));

    for (unsigned i = 0; i < length(pSiteFeatures.mEffectiveSites); ++i) {

        if (!pSiteFeatures.mEffectiveSites[i]) {
            continue;
        }

        k = 0;
        for (unsigned j = 0; j < length(seedTypes[i]); ++j) {
            mInputVec[k] = seedTypes[i][j];
            ++k;
        }
        for (unsigned j = 0; j < length(similarities[i]); ++j) {
            mInputVec[k] = similarities[i][j];
            ++k;
        }
        for (unsigned j = 0; j < length(auRichUp[i]); ++j) {
            mInputVec[k] = auRichUp[i][j];
            ++k;
        }
        for (unsigned j = 0; j < length(auRichDown[i]); ++j) {
            mInputVec[k] = auRichDown[i][j];
            ++k;
        }
        mInputVec[k] = sitePos[i][0];
        ++k;
        for (unsigned j = 0; j < length(seqMatch[i]); ++j) {
            mInputVec[k] = seqMatch[i][j];
            ++k;
        }
        mInputVec[k] = a1Match[i][0];

        mInputVec.normalize();
        mScores[i] = mModel.calc_score(mInputVec);

    }

    return 0;
}

template<class TRNAString>
void TSSVMSiteInputVector<TRNAString>::print_input_vector() {
    std::cout << mInputVec.transpose();
    std::cout << std::endl;
}

// Explicit template instantiation
template
class TSSVMSiteModel<TRNATYPE>;

template
class TSSVMSiteInputVector<TRNATYPE>;

} // namespace tssvm
