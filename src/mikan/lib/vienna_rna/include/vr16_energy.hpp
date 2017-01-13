#ifndef VR16_ENERGY_HPP_
#define VR16_ENERGY_HPP_

#include <string>
#include <map>

namespace vr16 {

//
// Energy parameters
//
class VR16EnergyParams
{
public:
    // Constant values
    const double K0;
    const double GASCONST;                     /* in [cal/K] */
    const int MAX_NINIO;                       /* maximum correction */
    static const int INF = 1000000;
    static const int NBPAIRS = 7;
    static const int NST = 0;                  /* Energy for nonstandard stacked pairs */
    static const int DEF = -50;                /* Default terminal mismatch, used for I  and any non_pairing bases */
    static const int NSM = 0;                  /* terminal mismatch for non standard pairs */

    // Declare variables
    double mTmeasure;                          /* temperature of param measurements */
    double mlxc37;                             /* parameter for logarithmic loop  energy extrapolation  */

    int **mStack37;
    int **mEnthalpies;                         /* stack enthalpies */

    int *mHairpin37;
    int *mBulge37;
    int *mInternalLoop37;

    int ***mMismatchI37;                       /* interior loop mismatches */
    int ***mMismatchH37;                       /* same for hairpins */
    int ***mMismatchM37;                       /* same for multiloops */
    int ***mMismH;                             /* mismatch enthalpies */

    int **mDangle5_37;                         /* 5' dangle exterior of pair */
    int **mDangle3_37;                         /* 3' dangle */
    int **mDangle3_H;                          /* corresponding enthalpies */
    int **mDangle5_H;

    /* constants for linearly destabilizing contributions for multi-loops
       F = ML_closing + ML_intern*(k-1) + ML_BASE*u  */
    int mMLBase37;
    int mMLClosing37;
    int mMLIntern37;

    /* Ninio-correction for asymmetric internal loops with branches n1 and n2 */
    /*    ninio_energy = min{max_ninio, |n1-n2|*F_ninio[min{4.0, n1, n2}] } */
    int mFNinio37[5];

    /* penalty for helices terminated by AU (actually not GC) */
    int mTerminalAU;

    /* penalty for forming bi-molecular duplex */
    int mDuplexInit;

    /* stabilizing contribution due to special hairpins of size 4 (tetraloops) */
    std::map<std::string, int> mTetraloops;    /* string containing the special tetraloops */
    int *mTetraEnergy37;                       /* Bonus energy for special tetraloops */
    int mTetraEnth37;

    std::map<std::string, int> mTriloops;      /* string containing the special triloops */
    int *mTriloopE37;                          /* Bonus energy for special Triloops */

    int ****mInt11_37;                         /* 1x1 interior loops */
    int ****mInt11_H;

    int *****mInt21_37;                        /* 2x1 interior loops */
    int *****mInt21_H;

    int ******mInt22_37;                       /* 2x2 interior loops */
    int ******mInt22_H;

public:
    // Define methods
    VR16EnergyParams(): K0(273.15), GASCONST(1.98717), MAX_NINIO(300)
    {
        init_heap();
        init_parameters();
    }
    ~VR16EnergyParams()
    {
        free_heap();
    }

private:
    void init_heap();
    void free_heap();

    void init_parameters();

    void init_stack37();
    void init_enthalpies();

    void init_hairpin37();
    void init_bulge37();
    void init_internal_loop37();

    void init_mismatch_i37();
    void init_mismatch_h37();
    void init_mismatch_m37();
    void init_mism_h();

    void init_dangle5_37();
    void init_dangle3_37();
    void init_dangle3_h();
    void init_dangle5_h();

    void init_tetraloops();
    void init_tetra_energy37();
    void init_triloop_e37();

    void init_int11_37();
    void init_int11_h();

    void init_int21_37();
    void init_int21_h();

    void init_int22_37();
    void init_int22_h();
};

} // namespace vr16

#endif /* VR16_ENERGY_HPP_ */
