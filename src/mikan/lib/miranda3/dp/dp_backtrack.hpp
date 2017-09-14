#ifndef DP_BACKTRACK_HPP_
#define DP_BACKTRACK_HPP_

#include <string>
#include <deque>
#include <algorithm>    // std::reverse
#include "dp_table.hpp"
#include "dp_score.hpp"
#include "dp_align.hpp"

namespace mr3dp {

class MR3DPTraceBack {
public:

    // Constructor
    MR3DPTraceBack() {};

    // Start backtracking
    void run(MR3DPTab &pTab, MR3DPScore &pScore, std::string &pQSeq, std::string &pDSeq) {
        // Clear path
        clear_path();

        int maxI = pTab.get_max_i();
        int maxJ = pTab.get_max_j();
        char maxTab = pTab.get_max_tab();

        mTerminate = false;
        if (maxTab == 'E') {
            back_track_e(pTab, pScore, pQSeq, pDSeq, maxI, maxJ);
        } else if (maxTab =='F') {
            back_track_f(pTab, pScore, pQSeq, pDSeq, maxI, maxJ);
        } else if (maxTab == 'G') {
            back_track_g(pTab, pScore, pQSeq, pDSeq, maxI, maxJ);
        }
    }

    void print_align() {
        mAlign.print_align();
    }

    // Get query sequence of alignment
    std::string &get_q_align() {
        return mAlign.get_q_align();
    }

    // Get database sequence of alignment
    std::string &get_d_align() {
        return mAlign.get_d_align();
    }

    int get_gap_q_count() {
        return mAlign.get_gap_q_count();
    }

    int get_gap_d_count() {
        return mAlign.get_gap_d_count();
    }

    void pirnt_path() {
        std::cout << "PathI: ";
        for (unsigned k = 0; k < mPathI.size(); k++) {
            std::cout << mPathI[k] << ", ";
        }
        std::cout << std::endl;

        std::cout << "PathJ: ";
        for (unsigned k = 0; k < mPathJ.size(); k++) {
            std::cout << mPathJ[k] << ", ";
        }
        std::cout << std::endl;
    }

private:
    bool mTerminate;

    // Path
    std::deque<int> mPathI;
    std::deque<int> mPathJ;

    MR3DPAlign mAlign;

    // Clear path
    void clear_path() {
        mPathI.clear();
        mPathJ.clear();
    }

    // Add element to path
    void add_to_path(int pI, int pJ) {
        // Add cell to path
        mPathI.push_back(pI);
        mPathJ.push_back(pJ);
    }

    // Remove the last element from path
    void remove_from_path(int pN) {
        // Remove cell from path
        while (pN > 0) {
            mPathI.pop_back();
            mPathJ.pop_back();
            --pN;
        }
    }

    // Backtrack - E
    void back_track_e(MR3DPTab &pTab, MR3DPScore &pScore, std::string &pQSeq, std::string &pDSeq, int pI, int pJ) {
        // Scores
        int gapE;
        int gapOpen;

        // Add the current cell to path
        add_to_path(pI, pJ);

        // Stop backtracking
        if (pI == 0 || pJ == 0) {
            remove_from_path(1);
            return;
        }

        // Gap penalties
        if (pI == 1) {
            gapE = pTab.get_e_val(pI, pJ) + pScore.gap_open();
        } else {
            gapE = pTab.get_e_val(pI, pJ) + pScore.gap_extend();
        }
        gapOpen =  pTab.get_e_val(pI, pJ) + pScore.gap_open();

        // Backtrack - E
        if (gapE ==  pTab.get_e_val(pI - 1, pJ)) {
            back_track_e(pTab, pScore, pQSeq, pDSeq, pI - 1, pJ);
            if (mTerminate) {
                return;
            }
        }

        // Backtrack - F
        if (gapOpen ==  pTab.get_f_val(pI - 1, pJ)) {
            back_track_f(pTab, pScore, pQSeq, pDSeq, pI - 1, pJ);
            if (mTerminate) {
                return;
            }
        }

        // Backtrack - G
        if (gapOpen ==  pTab.get_g_val(pI - 1, pJ)) {
            back_track_g(pTab, pScore, pQSeq, pDSeq, pI - 1, pJ);
            if (mTerminate) {
                return;
            }
        }

        // Remove the current cell from path
        remove_from_path(1);
    }

    // Backtrack - F
    void back_track_f(MR3DPTab &pTab, MR3DPScore &pScore, std::string &pQSeq, std::string &pDSeq, int pI, int pJ) {
        // Scores
        int gapF;
        int gapOpen;

        // Add the current cell to path
        add_to_path(pI, pJ);

        // Stop backtracking
        if (pI == 0 || pJ == 0) {
            remove_from_path(1);
            return;
        }

        // Gap penalties
        if (pJ == 1) {
            gapF =  pTab.get_f_val(pI, pJ) + pScore.gap_open();
        } else {
            gapF =  pTab.get_f_val(pI, pJ) + pScore.gap_extend();
        }
        gapOpen =  pTab.get_f_val(pI, pJ) + pScore.gap_open();

        // Backtrack - E
        if (gapOpen ==  pTab.get_e_val(pI, pJ - 1)) {
            back_track_e(pTab, pScore, pQSeq, pDSeq, pI, pJ - 1);
            if (mTerminate) {
                return;
            }
        }

        // Backtrack - F
        if (gapF ==  pTab.get_f_val(pI, pJ - 1)) {
            back_track_f(pTab, pScore, pQSeq, pDSeq, pI, pJ - 1);
            if (mTerminate) {
                return;
            }
        }

        // Backtrack - G
        if (gapOpen ==  pTab.get_g_val(pI, pJ - 1)) {
            back_track_g(pTab, pScore, pQSeq, pDSeq, pI, pJ - 1);
            if (mTerminate) {
                return;
            }
        }

        // Remove the current cell from path
        remove_from_path(1);
    }

    // Backtrack - G
    void back_track_g(MR3DPTab &pTab, MR3DPScore &pScore, std::string &pQSeq, std::string &pDSeq, int pI, int pJ) {
        // Score
        int scoreDiag;

        // Add the current cell to path
        add_to_path(pI, pJ);

        // Extend path to cell(0, 0) and add alignment
        if (pI == 0 || pJ == 0) {
            stop_bt(pQSeq, pDSeq, pI, pJ);
            mTerminate = true;
            return;
        }

        // Score
        scoreDiag = pTab.get_g_val(pI, pJ) - pScore.score_ab(pQSeq[pI - 1], pDSeq[pJ - 1]);

        // Backtrack - E
        if (scoreDiag ==  pTab.get_e_val(pI - 1, pJ - 1)) {
            back_track_e(pTab, pScore, pQSeq, pDSeq, pI - 1, pJ - 1);
            if (mTerminate) {
                return;
            }
        }

        // Backtrack - F
        if (scoreDiag ==  pTab.get_f_val(pI - 1, pJ - 1)) {
            back_track_f(pTab, pScore, pQSeq, pDSeq, pI - 1, pJ - 1);
            if (mTerminate) {
                return;
            }
        }

        // Backtrack - G
        if (scoreDiag ==  pTab.get_g_val(pI - 1, pJ - 1)) {
            back_track_g(pTab, pScore, pQSeq, pDSeq, pI - 1, pJ - 1);
            if (mTerminate) {
                return;
            }
        }

        // Remove the current cell from path
        remove_from_path(1);
    }

    // Stop backtracking
    void stop_bt(std::string &pQSeq, std::string &pDSeq, int pI, int pJ) {
        // Number of cells temporarily added to path
        int numCells;

        numCells = extend_path(pI, pJ);

        std::reverse(mPathI.begin(), mPathI.end());
        std::reverse(mPathJ.begin(), mPathJ.end());
        mAlign.create_align(pQSeq, pDSeq, mPathI, mPathJ);
        remove_from_path(numCells);
    }

    // Extend path to (0, 0)
    int extend_path(int pI, int pJ) {
        int n = 1;

        // Extend path to cell(0, 0)
        if (pI != 0) {
            for (int k = 0; k < pI; k++) {
                mPathI.push_back(pI - k - 1);
                mPathJ.push_back(0);
            }
            n += pI;
        } else if (pJ != 0) {
            for (int k = 0; k < pJ; k++) {
                mPathI.push_back(0);
                mPathJ.push_back(pJ - k - 1);
            }
            n += pJ;
        }

        return n;
    }

};

} // namespace mr3dp

#endif /* DP_BACKTRACK_HPP_ */