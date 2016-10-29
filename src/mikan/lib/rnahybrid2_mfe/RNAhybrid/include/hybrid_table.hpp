#ifndef HYBRID_TABLE_HPP_
#define HYBRID_TABLE_HPP_

#include <hybrid_energy.hpp>
#include <string>

namespace rh2{

const char TX = 5;

class RH2Table
{
public:
    void table_alloc(int max_target_len, int max_query_len);

    RH2Table(RH2EnergyFunc &pEn, int max_target_len, int max_query_len, std::string &p_seed_def):
        en(pEn),
        target_seq(0),
        query_seq(0),
        seed_def(p_seed_def),
        iloop_upper_limit(15),
        bloop_upper_limit(15),
        mTargetHelixStart(0),
        mTargetHelixEnd(0)
    {
        table_alloc(max_target_len, max_query_len);
    }

    int inpx(int i){return (int)(*target_seq)[i];}
    int inpy(int i){return (int)(*query_seq)[i];}

    void set_iloop_upper_limit(int l){iloop_upper_limit = l;}
    void set_bloop_upper_limit(int l){bloop_upper_limit = l;}

    void set_target_seq(std::vector<char> *seq, int qTargetLen);
    void set_target_seed_pos(int qTargetLen);
    void set_query_seq(std::vector<char> *seq, int qQueryLen);
    void set_query_seed_pos(int qQueryLen);

    int get_query_helix_start(){return mQueryHelixStart;}
    int get_query_helix_end(){return mQueryHelixEnd;}
    int get_target_helix_start(){return mTargetHelixStart;}
    int get_target_helix_end(){return mTargetHelixEnd;}

    float calc_hybrid(int i1, int j1, int i2, int j2);
    void calc_unpaired_left_top(int i1, int j1, int i2, int j2);
    void calc_unpaired_left_bot(int i1, int j1, int i2, int j2);
    void calc_closed(int i1, int j1, int i2, int j2);

    float get_unpaired_left_top(int i, int j){return tbl_unpaired_left_top[i][j];}
    float get_unpaired_left_bot(int i, int j){return tbl_unpaired_left_bot[i][j];}
    float get_closed(int i, int j){return tbl_closed[i][j];}

private:
    RH2EnergyFunc &en;

    std::vector<std::vector<float> > tbl_unpaired_left_top;
    std::vector<std::vector<float> > tbl_unpaired_left_bot;
    std::vector<std::vector<float> > tbl_closed;

    std::vector<char> *target_seq;
    std::vector<char> *query_seq;

    std::string &seed_def;

    int iloop_upper_limit, bloop_upper_limit;
    int seed_start, seed_end;
    int mQueryHelixStart, mQueryHelixEnd, mTargetHelixStart, mTargetHelixEnd;

};

} // namespace rh2

#endif /* HYBRID_TABLE_HPP_ */
