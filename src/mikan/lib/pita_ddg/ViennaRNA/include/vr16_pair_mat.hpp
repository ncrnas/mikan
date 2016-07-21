#ifndef VR16_PAIR_MAT_HPP_
#define VR16_PAIR_MAT_HPP_

namespace vr16 {

//
// Base pair matrix
//
class VR16PairMat
{
public:
    // Constant values
    static const int NBASES = 8;
    static const int MAXALPHA = 20;   /* maximal length of alphabet */

    // Declare variables
    std::string mLawAndOrder;
    int mBPPair[NBASES][NBASES];

    int mAlias[MAXALPHA+1];
    int mPair[MAXALPHA+1][MAXALPHA+1];
    int mRtype[8];

public:
    // Define methods
    VR16PairMat()
    {
        init_dat();
    }

    // Method prototypes
    int make_pair_matrix(int pEnergySet, std::string& pNonStandards, bool pNoGU);
    int encode_char(char pChr, int pEnergySet);
    int get_alias(int pIdx) {return mAlias[pIdx];}

private:
    void init_dat();

};

} // namespace vr16

#endif /* VR16_PAIR_MAT_HPP_ */
