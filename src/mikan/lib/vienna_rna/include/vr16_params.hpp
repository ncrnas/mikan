#ifndef VR16_PARAMS_HPP_
#define VR16_PARAMS_HPP_

#include <vr16_energy.hpp>                // VR16EnergyParams
#include <vr16_fold_options.hpp>          // VR16FoldOptions
#include <vr16_fold_options.hpp>          // VR16FoldOptions
#include <vr16_pair_mat.hpp>              // VR16PairMat
#include <cmath>
#include <map>

namespace vr16 {

//
// Scaled parameters for ViennaRNA fold & duplexfold
//
class VR16Params
{
public:
    // Constant values
    static const int NBPAIRS = 7;
    static const int MAXLOOP = 30; // must be MAXLOOP >= 30
    const double K0;
    const double GASCONST;                      /* in [cal/K] */
    const int INF;
    const int MAX_NINIO;                        /* maximum correction */

    // Declare variables
    int mId;
    int mStack[NBPAIRS+1][NBPAIRS+1];
    int mHairpin[31];
    int mBulge[MAXLOOP+1];
    int mInternalLoop[MAXLOOP+1];
    int mMismatchI[NBPAIRS+1][5][5];
    int mMismatchH[NBPAIRS+1][5][5];
    int mMismatchM[NBPAIRS+1][5][5];
    int mDangle5[NBPAIRS+1][5];
    int mDangle3[NBPAIRS+1][5];
    int mInt11[NBPAIRS+1][NBPAIRS+1][5][5];
    int mInt21[NBPAIRS+1][NBPAIRS+1][5][5][5];
    int mInt22[NBPAIRS+1][NBPAIRS+1][5][5][5][5];
    int mFNinio[5];
    double mlxc;
    int mMLBase;
    int mMLintern[NBPAIRS+1];
    int mMLclosing;
    int mTerminalAU;
    int mDuplexInit;
    int mTetraEnergy[200];
//    std::map<std::string, int> mTetraloops;
    int mTriloopE[40];
//    std::map<std::string, int> mTriloops;
    double mTemperature;

public:
    // Define methods
    VR16Params(VR16EnergyParams& pEn, double pTemprature) :
        K0(pEn.K0), GASCONST(pEn.GASCONST), INF(pEn.INF), MAX_NINIO(pEn.MAX_NINIO), mId(-1), mEn(pEn)
    {
        init_parameters(pTemprature);
    }

private:
    VR16EnergyParams& mEn;

private:
    void init_parameters(double pTemprature);

};

//
// Cache minimum values of internal loops
//
class VR16ParamIL
{

public:
    // Constant values
    static const int MAXLOOP = 30;
    static const int NBPAIRS = 7;

    // Define methods
    VR16ParamIL(VR16EnergyParams& pEn, VR16Params& pParams, VR16PairMat& pPairMat, VR16FoldOptions& pOpts) :
        mEn(pEn), mParams(pParams), mPairMat(pPairMat), mOpts(pOpts), INF(pParams.INF)
    {
        init_parameters();
    }
    int loop_energy(int pN1, int pN2, int pType1, int pType2, int pSi1, int pSj1, int pSp1, int pSq1);

public:
    int mMinIL[MAXLOOP+1][MAXLOOP+1];

private:
    VR16EnergyParams& mEn;
    VR16Params& mParams;
    VR16PairMat& mPairMat;
    VR16FoldOptions& mOpts;
    const int INF;
    static const int TURN = 3;

private:
    void init_parameters();
};

//
// Parameters of sub-optimal structure calculation
//
class VR16PFParams
{
public:
    // Constant values
    static const int NBPAIRS = 7;
    static const int MAXLOOP = 30;
    static const int PF_SCALE_10 = 10;
    static const int INITLEN = 50;

    // Declare variables
    int mInitLength;                   /* length in last call to init_pf_fold() */

    std::vector<double> mScale;

    std::vector<double> mExpHairpin;
    std::vector<double> mExpMLbase;

    double mExpBulge[MAXLOOP+1];
    double mExpInternal[MAXLOOP+1];

    double mExpNinio[5][MAXLOOP+1];

    double mExpMLclosing;
    double mExpMLintern[NBPAIRS+1];

    double mExpTermAU;
    double mExpDangle5[NBPAIRS+1][5];
    double mExpDangle3[NBPAIRS+1][5];

    double mExpTetra[40];
    double mExpTriloop[40];
    double mExpStack[NBPAIRS+1][NBPAIRS+1];

    double mExpMismatchI[NBPAIRS+1][5][5];
    double mExpMismatchH[NBPAIRS+1][5][5];
    double mExpMismatchM[NBPAIRS+1][5][5];

    double mExpInt11[NBPAIRS+1][NBPAIRS+1][5][5];
    double mExpInt21[NBPAIRS+1][NBPAIRS+1][5][5][5];
    double mExpInt22[NBPAIRS+1][NBPAIRS+1][5][5][5][5];

    double mExpMLContrib[NBPAIRS+1][5][5];

public:
    // Define methods
    VR16PFParams(VR16EnergyParams& pEn, double pTemperature, VR16Params& pParams, VR16FoldOptions& pOpts):
        mInitLength(0), mExpMLclosing(0), mExpTermAU(0), mInitTemp(0), mEn(pEn), mParams(pParams), mOpts(pOpts)
    {
        scale_pf_params(INITLEN, pTemperature);
    }

    // Method prototypes
    void scale_pf_params(unsigned int pLen, double pTemperature);
    void reset_scale();

private:
    double mInitTemp;                  /* temperature in last call to scale_pf_params */

    VR16EnergyParams& mEn;
    VR16Params& mParams;
    VR16FoldOptions& mOpts;

private:
    double pf_smooth(double pX);
};

} // namespace vr16

#endif /* VR16_PARAMS_HPP_ */
