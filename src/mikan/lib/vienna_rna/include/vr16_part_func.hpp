#ifndef VR16_PART_FUNC_HPP_
#define VR16_PART_FUNC_HPP_

#include <vr16_energy.hpp>                // VR16EnergyParams
#include <vr16_fold_options.hpp>          // VR16FoldOptions
#include <vr16_pair_mat.hpp>              // VR16PairMat
#include <vr16_params.hpp>                // VR16Params
#include <vector>
#include <string>

namespace vr16 {

//
// Calculate sub-optimal structures
//
class VR16PartFunc {
public:
    // Constant values
    static const int MAXLOOP = 30;
    static const int TURN = 3;
    static const int PLEN = 3;

    // Declare variables
    bool mStBack;

public:
    // Define methods
    VR16PartFunc(VR16EnergyParams &pEn, VR16Params &pParams, VR16PFParams &pPFParams, VR16PairMat &pPairMat,
                 VR16FoldOptions &pOpts) :
            mStBack(false), mInitLength(0), mEn(pEn), mParams(pParams), mPFParams(pPFParams), mPairMat(pPairMat),
            mOpts(pOpts) {}

    // Method prototypes
    float pf_fold(std::string &pString, std::string &pStructure, double pTemperature);

    void init_pf_fold(int pLen, double pTemperature);

private:
    int mInitLength;                   /* length in last call to init_pf_fold() */

    std::vector<int> mS;
    std::vector<int> mS1;

    std::vector<double> mQ;
    std::vector<double> mQb;
    std::vector<double> mQm;
    std::vector<double> mQm1;
    std::vector<double> mQqm;
    std::vector<double> mQqm1;
    std::vector<double> mQq;
    std::vector<double> mQq1;

    std::vector<int> mPtype;           /* precomputed array of pair types */

    std::vector<double> mPrml;
    std::vector<double> mPrmL0;
    std::vector<double> mPrmL1;
    std::vector<double> mQ1k;
    std::vector<double> mQln;

    std::vector<int> mJindx;

    VR16EnergyParams &mEn;
    VR16Params &mParams;
    VR16PFParams &mPFParams;
    VR16PairMat &mPairMat;
    VR16FoldOptions &mOpts;

private:
    void resize_arrays(unsigned int pLen);

    void make_ptypes(std::string &pStructure);

    double exp_hairpin_energy(int pU, int pType, int pSi1, int pSj1, int pIdx, std::string &pString);

    double exp_loop_energy(int pU1, int pU2, int pType, int pType2, int pSi1, int pSj1, int pSp1, int pSq1);

    void sprintf_bppm(int pLen, std::string &pStructure);

    char bppm_symbol(std::vector<float> &pP);

    double get_qqx_val(int pQQId, int pIdx);

    double get_qqmx_val(int pQQMId, int pIdx);

    double get_prmlx_val(int pPrmLId, int pIdx);

    void set_qqx_val(int pQQId, int pIdx, double pNewVal);

    void set_qqmx_val(int pQQMId, int pIdx, double pNewVal);

    void set_prmlx_val(int pPrmLId, int pIdx, double pNewVal);

    void backtrack(int pStrLen, int pIdPrmL0, int pIdPrmL1);

};

} // namespace vr16

#endif /* VR16_PART_FUNC_HPP_ */
