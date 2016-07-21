#ifndef TSSVM_MRNA_SVM_HPP_
#define TSSVM_MRNA_SVM_HPP_

#include <mikan/lib/two_step_svm/include/tssvm_mrna_feature.hpp>    // TSSVMRNARawFeatures
#include <Eigen/Dense>
#include <seqan/sequence.h>


namespace tssvm{

//
// RNA level SVM model
//
template <class TRNAString>
class TSSVMRNAModel
{
public:
    // Define methods
    TSSVMRNAModel(): mB(-1.92467490317), mHyperPlane(35)
    {
        init_hyperplane();
    }

    // Method prototypes
    int init_model();
    float calc_score(Eigen::VectorXf& pInput);

private:
    seqan::CharString mModelPath;
    const float mB;
    Eigen::VectorXf mHyperPlane;

private:
    void init_hyperplane();

};

//
// Input vector for RNA level SVM
//
template <class TRNAString>
class TSSVMRNAInputVector
{
public:
    // Define types
    typedef seqan::StringSet<seqan::String<float> > TFeatSet;

public:
    // Define methods
    TSSVMRNAInputVector(): mModel(), mInputVec(35) {}
    const seqan::String<float>& get_scores(){return mScores;}

    // Method prototypes
    void clear_scores();
    int classify(TSSVMRNARawFeatures<TRNAString> &pRNAFeatures);

private:
    TSSVMRNAModel<TRNAString> mModel;
    Eigen::VectorXf mInputVec;
    seqan::String<float> mScores;

private:
    int calc_score(TSSVMRNARawFeatures<TRNAString> &pRNAFeatures);
    void print_input_vector();

};

} // namespace tssvm

#endif /* TSSVM_MRNA_SVM_HPP_ */
