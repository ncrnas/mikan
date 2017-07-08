#ifndef TSSVM_MRNA_SVM_HPP_
#define TSSVM_MRNA_SVM_HPP_

#include <Eigen/Dense>
#include <seqan/sequence.h>
#include "mk_typedef.hpp"            // TRNATYPE, TCharSet, TRNASet, TIndexQGram, TFinder
#include "tssvm_mrna_feature.hpp"    // TSSVMRNARawFeatures

namespace tssvm {

//
// RNA level SVM model
//
class TSSVMRNAModel {
public:
    // Define methods
    TSSVMRNAModel() : mB(-1.92467490317f), mHyperPlane(35) {
        init_hyperplane();
    }

    // Method prototypes
    int init_model();

    float calc_score(Eigen::VectorXf &pInput);

private:
    const float mB;
    Eigen::VectorXf mHyperPlane;

private:
    void init_hyperplane();

};

//
// Input vector for RNA level SVM
//
class TSSVMRNAInputVector {
public:
    // Define types
    typedef seqan::StringSet<seqan::String<float> > TFeatSet;

public:
    // Define methods
    TSSVMRNAInputVector() : mModel(), mInputVec(35) {}

    const seqan::String<float> &get_scores() { return mScores; }

    // Method prototypes
    void clear_scores();

    int classify(TSSVMRNARawFeatures &pRNAFeatures);

private:
    TSSVMRNAModel mModel;
    Eigen::VectorXf mInputVec;
    seqan::String<float> mScores;

private:
    int calc_score(TSSVMRNARawFeatures &pRNAFeatures);

    void print_input_vector();

};

} // namespace tssvm

#endif /* TSSVM_MRNA_SVM_HPP_ */
