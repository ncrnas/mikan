#ifndef VR16_PAIR_MAT_HPP_
#define VR16_PAIR_MAT_HPP_

#include <string>

namespace vr16 {

//
// Base pair matrix
//
class VR16PairMat {
public:
    // Constant values
    static const int NBASES = 8;
    static const int MAXALPHA = 20;   /* maximal length of alphabet */

    // Declare variables
    std::string mLawAndOrder;
    int **mBPPair;

    int *mAlias;
    int **mPair;
    int *mRtype;

public:
    // Define methods
    VR16PairMat() {
        init_heap();
        init_dat();
    }

    ~VR16PairMat() {
        free_heap();
    }

    // Method prototypes
    int make_pair_matrix(int pEnergySet, std::string &pNonStandards, bool pNoGU);

    int encode_char(char pChr, int pEnergySet);

    int get_alias(int pIdx) { return mAlias[pIdx]; }

private:
    void init_heap();

    void free_heap();

    void init_dat();

};

} // namespace vr16

#endif /* VR16_PAIR_MAT_HPP_ */
