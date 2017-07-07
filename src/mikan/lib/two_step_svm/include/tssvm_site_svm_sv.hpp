#ifndef TSSVM_SITE_SVM_SV_HPP_
#define TSSVM_SITE_SVM_SV_HPP_

#include <Eigen/Dense>
#include <vector>

namespace tssvm {

//
// Support vectors of site level SVM model
//
class TSSVMSiteModelSV {
public:
    // Constant values
    static const int INPUT_FEAT_NUM = 95;

public:
    // Define methods
    TSSVMSiteModelSV(Eigen::MatrixXf &pSVs) : mSVs(pSVs), mSVDefVals(15391), mZeroIndices(95, false) {
        init_zero_indices();
        set_sv_def_vals();
    }

    // Define methods
    int init_sv_matix();

    void print_zero_indices();

    void print_def_vals();

    void print_sv();

private:
    Eigen::MatrixXf &mSVs;
    Eigen::VectorXf mSVDefVals;
    std::vector<bool> mZeroIndices;

private:
    void init_zero_indices();

    void set_sv_def_vals();
};


} // namespace tssvm

#endif /* TSSVM_SITE_SVM_SV_HPP_ */
