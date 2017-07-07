#ifndef ENERGY_HPP_
#define ENERGY_HPP_

#include <vector>
#include <algorithm>
#include <cmath>
#include <iostream>

namespace rh2 {

class RH2EnergyArray {
public:
    RH2EnergyArray() {
        init_heap();
    }

    ~RH2EnergyArray() {
        free_heap();
    }

    static const int ALPHASIZE = 5;
    static const int A = 0;
    static const int C = 1;
    static const int G = 2;
    static const int U = 3;
    static const int N = 4;
    static const int X = 5;

    float **il_asym_ar;

    float ****stack_dg_ar;
    float ****tstackh_dg_ar;
    float ****tstacki_dg_ar;
    float ******hl_tetra_ar;

    float *hl_ent_ar;
    float *bl_ent_ar;
    float *il_ent_ar;

    float ***dr_dangle_dg_ar;
    float ***dl_dangle_dg_ar;

    float ******int11_ar;
    float *******int21_ar;
    float ********int22_ar;

    float ******il_do_ar;

    char canPair[ALPHASIZE + 1][ALPHASIZE + 1];

    void init_energies();

private:
    void init_heap();

    void free_heap();

    void init_il_asym_ar();

    void init_stack_dg_ar();

    void init_tstackh_dg_ar();

    void init_tstacki_dg_ar();

    void init_hl_tetra_ar();

    void init_hl_ent_ar();

    void init_bl_ent_ar();

    void init_il_ent_ar();

    void init_dr_dangle_dg_ar();

    void init_dl_dangle_dg_ar();

    void init_int11_ar();

    void init_int21_ar();

    void init_int22_ar();

    void init_canPair();

};

class RH2TempEnergyArray {
public:
    RH2TempEnergyArray() {
        init_heap();
    }

    ~RH2TempEnergyArray() {
        free_heap();
    }

    void clear_array();

    float ****tmp_do_id;
    bool ****flag_bit_do_id;

private:
    void init_heap();

    void free_heap();

};

class RH2EnergyFunc : RH2EnergyArray {
public:
    RH2EnergyFunc() :
            target_seq(0),
            query_seq(0),
            e(2.718281828459f),
            t(273.15f),
            temp(273.15f + 37.0f),
            r(8.3143f),
            wkn(0.83f),
            npp(0.2f),
            pbp(0.1f * 0.83f),
            mloop_close(4.6f),
            free_base_penalty(0.4f),
            helix_penalty(0.1f),
            mF0(0.0f) {
        init_energies();
    }

    void set_target_seq(std::vector<char> *seq) { target_seq = seq; }

    void set_query_seq(std::vector<char> *seq) { query_seq = seq; }

    int inpx(int i) { return (int) (*target_seq)[i]; }

    int inpy(int i) { return (int) (*query_seq)[i]; }

    char is_pair(int i, int j);

    float sr_energy(int i, int j);

    float dr_energy(int i, int j);

    float dli_energy(int i, int j);

    float dl_energy(int i, int j);

    float dri_energy(int i, int j);

    //	float dangles (int i, int j, int i2, int j2, int k, int l, int k2,int l2);
    float sspenalty(int a);

    float log_interp(int size);

    float hl_ent(int size);

    float bl_ent(int size);

    float il_ent(int size);

    int lengthOf(int i, int j);

    float il_asym(int sl, int sr);

    float top_stack(int lb, int rb);

    float bot_stack(int lb, int rb);

    float bl_stacking(int t, int b, int i, int j);

    float il_stack_open(int i, int j);

    float il_stack_close(int i, int j);

    float int_special(int i, int j, int t, int b);

    float do_il_special(int i, int j, int k, int l, int u, int v, float e);

    float do_il(int i, int j, int k, int l, int u, int v, float e);

private:
    std::vector<char> *target_seq;
    std::vector<char> *query_seq;

    const float e;
    const float t;
    const float temp;
    const float r;

    const float wkn;

    const float npp;
    const float pbp;

    const float mloop_close;
    const float free_base_penalty;
    const float helix_penalty;

    float mF0;

};

inline char RH2EnergyFunc::is_pair(int i, int j) {
    return canPair[i][j];
}

inline float RH2EnergyFunc::sr_energy(int i, int j) {
    return (stack_dg_ar[inpx(i)][inpx(i + 1)][inpy(j + 1)][inpy(j)]);
}

inline float RH2EnergyFunc::dr_energy(int i, int j) {
    return (dr_dangle_dg_ar[inpx(i)][inpy(j)][inpy(j - 1)]);
}

inline float RH2EnergyFunc::dli_energy(int i, int j) {
    return (dr_dangle_dg_ar[inpy(j)][inpx(i)][inpx(i + 1)]);
}

inline float RH2EnergyFunc::dl_energy(int i, int j) {
    return (dl_dangle_dg_ar[inpx(i - 1)][inpx(i)][inpy(j)]);
}

inline float RH2EnergyFunc::dri_energy(int i, int j) {
    return (dl_dangle_dg_ar[inpy(j + 1)][inpy(j)][inpx(i)]);
}

//inline float RH2EnergyFunc::dangles( int, int j, int, int j2, int k, int, int k2, int)
//{
//	return ((dli_energy(j, k + 1) + dri_energy(j2, k2 + 1)) * wkn);
//}

inline float RH2EnergyFunc::sspenalty(int a) {
    return (npp * a);
}

inline float RH2EnergyFunc::log_interp(int size) {
    return (1.079 * log(((float) size) / 30.0));
}

inline float RH2EnergyFunc::hl_ent(int size) {
    return (size <= 30 ? hl_ent_ar[size] : hl_ent_ar[30] + log_interp(size));
}

inline float RH2EnergyFunc::bl_ent(int size) {
    return (size <= 30 ? bl_ent_ar[size] : bl_ent_ar[30] + log_interp(size));
}

inline float RH2EnergyFunc::il_ent(int size) {
    return ((size) <= 30 ? il_ent_ar[(size)] : il_ent_ar[30] + log_interp((size)));
}

inline int RH2EnergyFunc::lengthOf(int i, int j) {
    return (j - i);
}

inline float RH2EnergyFunc::il_asym(int sl, int sr) {
    return (std::min(3.0, abs(sl - sr) * 0.3));
}

inline float RH2EnergyFunc::il_stack_open(int i, int j) {
    return tstacki_dg_ar[inpx(i)][inpx((i) + 1)][inpy((j) + 1)][inpy(j)];
}

inline float RH2EnergyFunc::il_stack_close(int i, int j) {
    return tstacki_dg_ar[inpy(j)][inpy((j) - 1)][inpx((i) - 1)][inpx(i)];
}

inline float RH2EnergyFunc::bl_stacking(int t, int b, int i, int j) {
    if ((t == 0) && (b == 1)) {
        return (stack_dg_ar[inpx(i)][inpx(i + 1)][inpy(j + 2)][inpy(j)]);
    } else if ((t == 1) && (b == 0)) {
        return (stack_dg_ar[inpx(i)][inpx(i + 2)][inpy(j + 1)][inpy(j)]);
    } else {
        return (mF0);
    }
}

inline float RH2EnergyFunc::int_special(int i, int j, int t, int b) {
    if ((t == 1) && (b == 1)) {
        return (int11_ar[inpx(i)][inpy(j)][inpx(i + 1)][inpy(j + 1)][inpx(i + 2)][inpy(j + 2)]);
    } else if ((t == 1) && (b == 2)) {
        return (int21_ar[inpx(i)][inpy(j)][inpx(i + 1)][inpy(j + 1)][inpy(j + 2)][inpx(i + 2)][inpy(j + 3)]);
    } else if ((t == 2) && (b == 1)) {
        return (int21_ar[inpy(j + 2)][inpx(i + 3)][inpy(j + 1)][inpx(i + 2)][inpx(i + 1)][inpy(j)][inpx(i)]);
    } else if ((t == 2) && (b == 2)) {
        return (int22_ar[inpx(i)][inpy(j)][inpx(i + 1)][inpx(i + 2)][inpy(j + 1)]
        [inpy(j + 2)][inpx(i + 3)][inpy(j + 3)]);
    } else {
        return (mF0);
    }
}

inline float RH2EnergyFunc::do_il_special(int i, int j, int k, int l, int u, int v, float e) {
    return (e + int_special(i, j, l - k, v - u));
}

inline float RH2EnergyFunc::do_il(int, int, int k, int l, int u, int v, float e) {

    return (e + il_stack_close(l + 1, v + 1) + il_ent(l - k + v - u) + il_asym(l - k, v - u));

}

} // namespace rh2

#endif /* ENERGY_HPP_ */
