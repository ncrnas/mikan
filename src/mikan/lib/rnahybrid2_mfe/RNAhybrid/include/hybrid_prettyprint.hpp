#ifndef HYBRID_PRETTYPRINT_HPP_
#define HYBRID_PRETTYPRINT_HPP_

#include <hybrid_backtrace.hpp>
#include <string>
#include <algorithm>

namespace rh2 {

class RH2PrettyPrinter {
public:

    void alloc_string(int targetlength, int querylength) {
        t1.resize(2 * std::max(targetlength, querylength), 0);
        t2.resize(2 * std::max(targetlength, querylength), 0);
        t3.resize(2 * std::max(targetlength, querylength), 0);
        t4.resize(2 * std::max(targetlength, querylength), 0);
        r1 = t1.begin();
        r2 = t2.begin();
        r3 = t3.begin();
        r4 = t4.begin();
    }

    RH2PrettyPrinter(RH2BackTrace &pBt, int targetlength, int querylength) :
            bt(pBt), a1(0), a2(0), a3(0), a4(0), a5(0), a6(0), target_seq(0), query_seq(0),
            mCountToA1(0) {
        alloc_string(targetlength, querylength);
    }

    int inpx(int i) { return (int) (*target_seq)[i]; }

    int inpy(int i) { return (int) (*query_seq)[i]; }

    void set_target_seq(std::vector<char> *seq) { target_seq = seq; }

    void set_query_seq(std::vector<char> *seq) { query_seq = seq; };

    void sx(std::string::iterator &t, int i);

    void shift(std::string::iterator &t, int u, int v);

    void one_space(std::string::iterator &t);

    void ssx(std::string::iterator &t, int i, int j);

    void sy(std::string::iterator &t, int i);

    void ssy(std::string::iterator &t, int i, int j);

    void blanks(std::string::iterator &t, int i, int j);

    void blanks2(std::string::iterator &t, int i, int j);

    void asym1(std::string::iterator &t, int l, int r, int u, int v);

    void asym2(std::string::iterator &t, int l, int r, int u, int v);

    void pp_str_Hybrid();

    int get_count_to_A1() { return mCountToA1; }

    int get_hit_length();

    void reset() {
        std::fill(t1.begin(), t1.end(), 0);
        std::fill(t2.begin(), t2.end(), 0);
        std::fill(t3.begin(), t3.end(), 0);
        std::fill(t4.begin(), t4.end(), 0);

        r1 = t1.begin();
        r2 = t2.begin();
        r3 = t3.begin();
        r4 = t4.begin();

        mCountToA1 = 0;
    }

    std::string t1, t2, t3, t4;

private:
    RH2BackTrace &bt;

    std::string::iterator r1, r2, r3, r4;

    int a1, a2, a3, a4, a5, a6;

    std::vector<char> *target_seq;
    std::vector<char> *query_seq;

    int mCountToA1;

};

} // namespace rh2


#endif /* HYBRID_PRETTYPRINT_HPP_ */
