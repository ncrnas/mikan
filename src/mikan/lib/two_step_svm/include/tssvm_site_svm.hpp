#ifndef TSSVM_SITE_SVM_HPP_
#define TSSVM_SITE_SVM_HPP_

#include <tssvm_site_feature.hpp>   // TSSVMRawFeatures
#include <Eigen/Dense>
#include <seqan/sequence.h>

namespace tssvm{

//
// Site level SVM model
//
template <class TRNAString>
class TSSVMSiteModel
{
public:
    // Constant values
    static const int INPUT_FEAT_NUM = 95;
    static const int MODEL_M = 15391;
    static const int K_POLY_DEGREE = 5;

public:
    // Define methods
    TSSVMSiteModel(): mB(-1.26691213422f), mAlphas(15391), mSVs(15391, 95), mMatProd(15391)
    {}

    // Define methods
    int init_model(seqan::CharString& pModelPath);
    float calc_score(Eigen::VectorXf& pInput);

private:
    seqan::CharString mModelPath;
    const float mB;
    Eigen::VectorXf mAlphas;
    Eigen::MatrixXf mSVs;
    Eigen::VectorXf mMatProd;

private:
    int init_alpha();
    int init_sv(seqan::CharString& pModelPath);

};

//
// Input vector for Site level SVM
//
template <class TRNAString>
class TSSVMSiteInputVector
{
public:
    // Define types
    typedef seqan::StringSet<seqan::String<float> > TFeatSet;

public:
    // Define methods
    TSSVMSiteInputVector(TSSVMSiteModel<TRNAString>& pModel): mModel(pModel), mInputVec(95) {}
    const seqan::String<float>& get_scores(){return mScores;}

    // Method prototypes
    void clear_scores();
    int classify(TSSVMRawFeatures<TRNAString> &pSiteFeatures);

private:
    TSSVMSiteModel<TRNAString>& mModel;
    Eigen::VectorXf mInputVec;
    seqan::String<float> mScores;

private:
    int calc_score(TSSVMRawFeatures<TRNAString> &pSiteFeatures);
    void print_input_vector();

};

} // namespace tssvm

#endif /* TSSVM_SITE_SVM_HPP_ */
