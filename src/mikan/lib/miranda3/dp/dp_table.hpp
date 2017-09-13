#ifndef DP_TABLE_HPP_
#define DP_TABLE_HPP_

#include <string>
#include <iostream>
#include <iomanip>      // std::setw
#include "mk_nd_array.hpp"

namespace mr3dp {

class MR3DPTab {

public:

    // Constructor
    MR3DPTab() {
        mTabN = 0;
        mTabM = 0;
        mSeqQN = 0;
        mSeqDN = 0;

        clear_max_vals();
    }

    ~MR3DPTab() {
        clear_tables();
    }

    void clear_tables() {
        if (mTabN != 0 && mTabM != 0) {
            mikan::delete_2d_array<int>(mTabE, mTabN);
            mikan::delete_2d_array<int>(mTabF, mTabN);
            mikan::delete_2d_array<int>(mTabG, mTabN);
        }

        mTabN = 0;
        mTabM = 0;
    }

    void reset_tables() {
        if (mTabN < mSeqQN || mTabM < mSeqDN) {
            clear_tables();

            mTabE = mikan::create_2d_array<int>(mSeqQN, mSeqDN);
            mTabF = mikan::create_2d_array<int>(mSeqQN, mSeqDN);
            mTabG = mikan::create_2d_array<int>(mSeqQN, mSeqDN);

            mTabN = mSeqQN;
            mTabM = mSeqDN;
        }
    }

    // Create DP table
    void create_tables(std::string &qQSeq, std::string& qDSeq) {
        mSeqQN =  qQSeq.size() + 1;
        mSeqDN = qDSeq.size() + 1;

        reset_tables();
        clear_max_vals();
    }

    // Print tables
    void print_tables(std::string &qQSeq, std::string& qDSeq) {
        std::cout << std::endl << "### Table E ###" << std::endl;
        print_tab(qQSeq, qDSeq, 'E');

        std::cout << "### Table F ###" << std::endl;
        print_tab(qQSeq, qDSeq, 'F');

        std::cout << "### Table G ###" << std::endl;
        print_tab(qQSeq, qDSeq, 'G');
    }

    // Row size
    int get_row_size() {
        return mTabN;
    }

    // Column size
    int get_col_size() {
        return mTabM;
    }

    // Update cell value
    void set_e_cell(int pI, int pJ, int pD) {
        set_cell(pI, pJ, pD, 'E');
    }

    void set_f_cell(int pI, int pJ, int pD) {
        set_cell(pI, pJ, pD, 'F');
    }

    void set_g_cell(int pI, int pJ, int pD) {
        set_cell(pI, pJ, pD, 'G');
    }

    // Get cell value
    int get_e_val(int i, int j) {
        return mTabE[i][j];
    }

    int get_f_val(int i, int j) {
        return mTabF[i][j];
    }

    int get_g_val(int i, int j) {
        return mTabG[i][j];
    }

    // Get max score
    int get_max_score() {
        return mMaxScore;
    }

    // Get i of the cell with max score
    int get_max_i() {
        return mMaxI;
    }

    // Get j of the cell with max score
    int get_max_j() {
        return mMaxJ;
    }

    // Get tab of max score
    char get_max_tab() {
        return mMaxTab;
    }


private:
    // Table
    int **mTabE;
    int **mTabF;
    int **mTabG;

    // Table size
    int mTabN;
    int mTabM;

    // Sequence size
    int mSeqQN;
    int mSeqDN;

    // Max score and cells
    char mMaxTab;
    int mMaxScore;
    int mMaxI;
    int mMaxJ;


    // Print table
    void print_tab(std::string &qQSeq, std::string& qDSeq, char pTab) {

        // Print database sequence
        std::cout << "         ";
        for (unsigned i = 0; i < qDSeq.size(); ++i) {
            std::cout << "      " << qDSeq[i];
        }
        std::cout << std::endl;

        // Print rows
        for (int i = 0; i < mTabN; ++i) {
            for (int j = 0; j < mTabM; ++j) {

                // Print query sequence
                if (j == 0) {
                    if (i == 0) {
                        std::cout << "  ";
                    } else {
                        std::cout << " " << qQSeq[i - 1];
                    }
                }

                // Print cell value
                if (pTab == 'E') {
                    std::cout  << std::setw(7) << mTabE[i][j];
                } else if (pTab == 'F') {
                    std::cout  << std::setw(7) << mTabF[i][j];
                } else if (pTab == 'G') {
                    std::cout  << std::setw(7) << mTabG[i][j];
                }

            }
            std::cout << std::endl;
        }
        std::cout << std::endl;
    }

    // Update cell value
    void set_cell(int i, int j, int d, char pTab) {
//        if (d < 0) {
//            d = 0;
//        }

        if (pTab == 'E') {
            mTabE[i][j] = d;
        } else if (pTab == 'F') {
            mTabF[i][j] = d;
        } else if (pTab == 'G') {
            mTabG[i][j] = d;
        }

        if (d >= mMaxScore) {
            mMaxScore = d;
            mMaxI = i;
            mMaxJ = j;
            mMaxTab = pTab;
        }
    }

    void clear_max_vals() {
        mMaxScore = 0;
        mMaxI = 0;
        mMaxJ = 0;
        mMaxTab = 'G';
    }
};

} // namespace mr3dp

#endif /* DP_TABLE_HPP_ */
