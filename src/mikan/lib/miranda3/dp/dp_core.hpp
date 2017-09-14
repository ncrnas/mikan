#ifndef DP_CORE_HPP_
#define DP_CORE_HPP_

#include <string>
#include "dp_table.hpp"
#include "dp_score.hpp"
#include "dp_backtrack.hpp"

namespace mr3dp {

class MR3DPCore {

public:
    // Constructor
    MR3DPCore() {}

    // Start table update
    void run(std::string &pQSeq, std::string &pDSeq, int pV00) {
        mTab.create_tables(pQSeq, pDSeq);
        init_table(pV00);
        update_cells(pQSeq, pDSeq);

        mBT.run(mTab, mScore, pQSeq, pDSeq);

//        mTab.print_tables(pQSeq, pDSeq);
//        mBT.print_align();
    }

    // Get max score
    int get_max_score() {
        return mTab.get_max_score();
    }

    // Get query sequence of alignment
    std::string &get_q_align() {
        return mBT.get_q_align();
    }

    // Get database sequence of alignment
    std::string &get_d_align() {
        return mBT.get_d_align();
    }

    int get_gap_q_count() {
        return mBT.get_gap_q_count();
    }

    int get_gap_d_count() {
        return mBT.get_gap_d_count();
    }

private:
    // DP table
    MR3DPTab mTab;

    // Scoring scheme
    MR3DPScore mScore;

    // Process backtracking
    MR3DPTraceBack mBT;

    // Initialize DP tables
    void init_table(int pV00) {
//        pV00 = 0;

        // Cell (0, 0)
        mTab.set_e_cell(0, 0, pV00);
        mTab.set_f_cell(0, 0, pV00);
        mTab.set_g_cell(0, 0, pV00);

        int gap;
        // First column
        for (int i = 1; i < mTab.get_col_size(); i++) {
            gap = mScore.gap_open() + (i - 1) * mScore.gap_extend();
            mTab.set_e_cell(0, i, pV00 - gap);
            mTab.set_f_cell(0, i, pV00 - gap);
            mTab.set_g_cell(0, i, pV00 - gap);
        }

        // First row
        for (int i = 1; i < mTab.get_row_size(); i++) {
            gap = mScore.gap_open() + (i - 1) * mScore.gap_extend();
            mTab.set_e_cell(i, 0, pV00 - gap);
            mTab.set_f_cell(i, 0, pV00 - gap);
            mTab.set_g_cell(i, 0, pV00 - gap);
        }
    }

    // Update cells
    void update_cells(std::string &pQSeq, std::string &pDSeq) {
        // Scores
        int scoreHoriz;
        int scoreVert;
        int scoreDiag;

        // Gap penalties
        int gapE;
        int gapF;

        // Score for Match/mismatch score
        int scoreAB;
        int cellScore;

        // Update cells
        for (int i = 1; i < mTab.get_row_size(); ++i) {
            for (int j = 1; j < mTab.get_col_size(); ++j) {
                // Match/mismatch score
                scoreAB = mScore.score_ab(pQSeq[i - 1], pDSeq[j - 1]);

                // Gap penalties for E and F
                if (i == 1) {
                    gapE = mScore.gap_open();
                } else {
                    gapE = mScore.gap_extend();
                }
                if (j == 1) {
                    gapF = mScore.gap_open();
                } else {
                    gapF = mScore.gap_extend();
                }

                // Scores for updating table E (vertical)
                scoreVert = mTab.get_e_val(i - 1, j) - gapE;
                scoreHoriz = mTab.get_f_val(i - 1, j) - mScore.gap_open();
                scoreDiag = mTab.get_g_val(i - 1, j) - mScore.gap_open();

                // Update E (tabE)
                cellScore = (scoreVert > scoreHoriz) ? scoreVert : scoreHoriz;
                cellScore = (scoreDiag > cellScore) ? scoreDiag : cellScore;
                mTab.set_e_cell(i, j, cellScore);

                // Scores for updating table F (horizontal)
                scoreVert = mTab.get_e_val(i, j - 1) - mScore.gap_open();
                scoreHoriz = mTab.get_f_val(i, j - 1) - gapF;
                scoreDiag = mTab.get_g_val(i, j - 1) - mScore.gap_open();

                // Update F (tabF)
                cellScore = (scoreVert > scoreHoriz) ? scoreVert : scoreHoriz;
                cellScore = (scoreDiag > cellScore) ? scoreDiag : cellScore;
                mTab.set_f_cell(i, j, cellScore);

                // Scores for updating table G (diagonal)
                scoreVert = mTab.get_e_val(i - 1, j - 1) + scoreAB;
                scoreHoriz = mTab.get_f_val(i - 1, j - 1) + scoreAB;
                scoreDiag = mTab.get_g_val(i - 1, j - 1) + scoreAB;

                // Update G (tabG)
                cellScore = (scoreVert > scoreHoriz) ? scoreVert : scoreHoriz;
                cellScore = (scoreDiag > cellScore) ? scoreDiag : cellScore;
                mTab.set_g_cell(i, j, cellScore);

            }
        }
    }
};

} // namespace mr3dp

#endif /* DP_CORE_HPP_ */