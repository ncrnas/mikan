#ifndef TM1_MRNA_SVM_HPP_
#define TM1_MRNA_SVM_HPP_

#include <Eigen/Dense>
#include <seqan/sequence.h>
#include "mk_typedef.hpp"         // TRNATYPE, TCharSet, TRNASet, TIndexQGram, TFinder
#include "tm1_mrna_feature.hpp"   // TM1ScaledFeatures,

namespace tm1p {

//
// mRNA level SVM model
//
class TM1MRNAModel {
public:
    // Constant values
    static const unsigned INPUT_FEAT_NUM = 30;
    static const unsigned MODEL_M = 435;

public:
    // Define methods
    TM1MRNAModel() : mRho(-0.279867f), mGamma(0.09f), mAlphas(435), mSVs(435, 30), mSquaredSVs(435), mMatProd(435) {
        init_model();
    }

    // Define methods
    float calc_score(Eigen::VectorXf &pInput, float pValSquared);

private:
    const float mRho;
    const float mGamma;
    Eigen::VectorXf mAlphas;
    Eigen::MatrixXf mSVs;
    Eigen::VectorXf mSquaredSVs;
    Eigen::VectorXf mMatProd;

private:
    void init_model();

    void init_alpha();

    void init_sv();

    void init_sv_squared();

};

//
// Input vector for mRNA level SVM
//
class TM1MRNAInputVector {
public:
    // Constant values
    static const unsigned INPUT_FEAT_NUM = 30;

public:
    // Define methods
    TM1MRNAInputVector() : mInputVec(30) {}

    const seqan::String<float> &get_scores() { return mScores; }

    // Method prototypes
    void clear_scores();

    int classify(TM1ScaledFeatures  &pMRNAFeatures);

private:
    TM1MRNAModel mModel;
    Eigen::VectorXf mInputVec;
    seqan::String<float> mScores;

private:
    int calc_score(TM1ScaledFeatures  &pMRNAFeatures);

    void print_input_vector(float pValSquared);

};

} // namespace tm1p

#endif /* TM1_MRNA_SVM_HPP_ */
