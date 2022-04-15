#include "gaussian_elimination.hpp"

Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> mod2Echelonize(Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> &m) {
    Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> linDependencies = Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>::Identity(m.rows(), m.rows());

    std::size_t currentNonReducedRow = 0;
    std::size_t currentPivotIndex;
    std::size_t nextPivotIndex;

    for (std::size_t col = 0; col < m.cols(); col++) {
        nextPivotIndex = std::numeric_limits<std::size_t>::max();

        // get the pivot row for column col
        for (std::size_t i = currentNonReducedRow; i < m.rows(); i++) {
            if (m(i, col)) {
                nextPivotIndex = i;
                break;
            }
        }

        // no candidate row found, continue to next column
        if (nextPivotIndex == std::numeric_limits<std::size_t>::max()) {
            continue;
        }

        // the pivot is different from the current non-reduced row, switch them
        if (nextPivotIndex != currentNonReducedRow) {
            m.row(nextPivotIndex).swap(m.row(currentNonReducedRow));
            linDependencies.row(nextPivotIndex).swap(linDependencies.row(currentNonReducedRow));
        }
        currentPivotIndex = currentNonReducedRow;
        currentNonReducedRow++;

        // reduce all the remaining rows with the pivot
        for (std::size_t i = currentNonReducedRow; i < m.rows(); i++) {
            if (m(i, col)) {
                // xor
                for (std::size_t j = 0; j < m.cols(); j++) {
                    m(i, j) = m(i, j) ^ m(currentPivotIndex, j);
                }
                for (std::size_t j = 0; j < linDependencies.cols(); j++) {
                    linDependencies(i, j) = linDependencies(i, j) ^ linDependencies(currentPivotIndex, j);
                }
            }
        }
    }

    std::size_t lastNonZeroRow = linDependencies.rows() - 1;
    std::size_t numLinRelations = 0;

    while (lastNonZeroRow >= 0 && m.row(lastNonZeroRow).isZero()) {
        lastNonZeroRow--;
        numLinRelations++;
    }

    return linDependencies(Eigen::seqN(lastNonZeroRow + 1, numLinRelations), Eigen::all);
}
