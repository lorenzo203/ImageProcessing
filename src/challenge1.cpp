#include <Eigen/Dense>
#include <Eigen/IterativeLinearSolvers>
#include <Eigen/Sparse>
#include <unsupported/Eigen/SparseExtra>
#include "lis.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

namespace {

using ImageMatrix = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
using ByteImageMatrix = Eigen::Matrix<unsigned char, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;

unsigned char clamp_pixel(const double value) {
    const double rounded = std::round(value);
    const double clamped = std::clamp(rounded, 0.0, 255.0);
    return static_cast<unsigned char>(clamped);
}

ByteImageMatrix vector_to_image(const Eigen::VectorXd &vector, const int rows, const int cols) {
    ByteImageMatrix image(rows, cols);
    for (int r = 0; r < rows; ++r) {
        for (int c = 0; c < cols; ++c) {
            const Eigen::Index idx = static_cast<Eigen::Index>(r) * static_cast<Eigen::Index>(cols) + c;
            image(r, c) = clamp_pixel(vector(idx));
        }
    }
    return image;
}

void save_png(const std::string &filename, const ByteImageMatrix &image) {
    if (!stbi_write_png(filename.c_str(), static_cast<int>(image.cols()), static_cast<int>(image.rows()), 1,
                        image.data(), static_cast<int>(image.cols()))) {
        throw std::runtime_error("Failed to write image: " + filename);
    }
}

void export_vector_market(const std::string &filename, const Eigen::VectorXd &vector) {
    std::ofstream out(filename);
    if (!out) {
        throw std::runtime_error("Unable to open " + filename + " for writing");
    }
    out << "%%MatrixMarket vector coordinate real general\n";
    out << vector.size() << '\n';
    for (Eigen::Index i = 0; i < vector.size(); ++i) {
        out << (i + 1) << ' ' << vector(i) << '\n';
    }
}

Eigen::MatrixXd create_filter_hav1() {
    Eigen::MatrixXd kernel(3, 3);
    kernel << 1.0, 1.0, 0.0,
              1.0, 2.0, 1.0,
              0.0, 1.0, 1.0;
    return kernel * (1.0 / 8.0);
}

Eigen::MatrixXd create_filter_hsh1() {
    Eigen::MatrixXd kernel(3, 3);
    kernel << 0.0, -2.0, 0.0,
             -2.0,  9.0, -2.0,
              0.0, -2.0, 0.0;
    return kernel;
}

Eigen::MatrixXd create_filter_hed2() {
    Eigen::MatrixXd kernel(3, 3);
    kernel << -1.0, -2.0, -1.0,
               0.0,  0.0,  0.0,
               1.0,  2.0,  1.0;
    return kernel;
}

Eigen::SparseMatrix<double> create_convolution_matrix(const Eigen::MatrixXd &kernel, const int rows, const int cols) {
    if (kernel.rows() != kernel.cols() || kernel.rows() % 2 == 0) {
        throw std::invalid_argument("Kernel must be square with odd dimension");
    }

    const int size = static_cast<int>(kernel.rows());
    const int radius = size / 2;
    const long total_size = static_cast<long>(rows) * cols;

    std::vector<Eigen::Triplet<double>> triplets;
    triplets.reserve(static_cast<size_t>(total_size) * static_cast<size_t>(size * size / 2));

    for (int r = 0; r < rows; ++r) {
        for (int c = 0; c < cols; ++c) {
            const long row_index = static_cast<long>(r) * cols + c;
            for (int i = -radius; i <= radius; ++i) {
                const int rr = r + i;
                if (rr < 0 || rr >= rows) {
                    continue;
                }
                for (int j = -radius; j <= radius; ++j) {
                    const int cc = c + j;
                    if (cc < 0 || cc >= cols) {
                        continue;
                    }
                    const double value = kernel(i + radius, j + radius);
                    if (value == 0.0) {
                        continue;
                    }
                    const long col_index = static_cast<long>(rr) * cols + cc;
                    triplets.emplace_back(row_index, col_index, value);
                }
            }
        }
    }

    Eigen::SparseMatrix<double> convolution(total_size, total_size);
    convolution.setFromTriplets(triplets.begin(), triplets.end());
    return convolution;
}

struct LisSolveResult {
    bool success{false};
    LIS_INT iterations{0};
    double residual{0.0};
};

LisSolveResult solve_with_lis(const Eigen::SparseMatrix<double> &matrix, const Eigen::VectorXd &rhs, Eigen::VectorXd &solution) {
    LisSolveResult result{};

    int argc = 1;
    char program_name[] = "challenge1";
    char *argv_values[] = {program_name, nullptr};
    char **argv = argv_values;
    lis_initialize(&argc, &argv);

    Eigen::SparseMatrix<double, Eigen::RowMajor> row_major = matrix;
    row_major.makeCompressed();

    const LIS_INT n = static_cast<LIS_INT>(row_major.rows());
    const LIS_INT nnz = static_cast<LIS_INT>(row_major.nonZeros());

    LIS_INT *ptr = nullptr;
    LIS_INT *index = nullptr;
    LIS_SCALAR *value = nullptr;

    const LIS_INT alloc_status = lis_matrix_malloc_csr(n, nnz, &ptr, &index, &value);
    if (alloc_status != LIS_SUCCESS) {
        lis_finalize();
        return result;
      }

    const int *outer = row_major.outerIndexPtr();
    const int *inner = row_major.innerIndexPtr();
    const double *vals = row_major.valuePtr();

    for (LIS_INT i = 0; i <= n; ++i) {
        ptr[i] = static_cast<LIS_INT>(outer[i]);
    }
    for (LIS_INT k = 0; k < nnz; ++k) {
        index[k] = static_cast<LIS_INT>(inner[k]);
        value[k] = static_cast<LIS_SCALAR>(vals[k]);
    }

    LIS_MATRIX A;
    lis_matrix_create(LIS_COMM_WORLD, &A);
    lis_matrix_set_size(A, 0, n);
    lis_matrix_set_type(A, LIS_MATRIX_CSR);
    lis_matrix_set_csr(nnz, ptr, index, value, A);
    lis_matrix_assemble(A);

    LIS_VECTOR b;
    LIS_VECTOR x;
    lis_vector_create(LIS_COMM_WORLD, &b);
    lis_vector_set_size(b, 0, n);
    lis_vector_duplicate(b, &x);

    for (LIS_INT i = 0; i < n; ++i) {
        lis_vector_set_value(LIS_INS_VALUE, i, static_cast<LIS_SCALAR>(rhs[static_cast<Eigen::Index>(i)]), b);
    }

    LIS_SOLVER solver;
    lis_solver_create(&solver);
    lis_solver_set_option(const_cast<char *>("-i cgs"), solver);
    lis_solver_set_option(const_cast<char *>("-p ilu"), solver);
    lis_solver_set_option(const_cast<char *>("-tol 1.0e-14"), solver);
    lis_solver_set_option(const_cast<char *>("-maxiter 5000"), solver);

    const LIS_INT status = lis_solve(A, b, x, solver);
    result.success = (status == LIS_SUCCESS);

    lis_solver_get_iter(solver, &result.iterations);
    lis_solver_get_residualnorm(solver, &result.residual);

    solution.resize(static_cast<Eigen::Index>(n));
    for (LIS_INT i = 0; i < n; ++i) {
        LIS_SCALAR value_i = 0.0;
        lis_vector_get_value(x, i, &value_i);
        solution[static_cast<Eigen::Index>(i)] = static_cast<double>(value_i);
    }

    lis_solver_destroy(solver);
    lis_vector_destroy(b);
    lis_vector_destroy(x);
    lis_matrix_destroy(A);
    lis_finalize();

    return result;
}

} // namespace

int main() {
    try {
        int width = 0;
        int height = 0;
        int channels = 0;
        unsigned char *image_data = stbi_load("data/uma.jpg", &width, &height, &channels, 1);
        if (!image_data) {
            throw std::runtime_error("Unable to load uma.jpg");
        }

        const Eigen::Map<const ByteImageMatrix> input_view(image_data, height, width);
        ByteImageMatrix original_uint = input_view;
        stbi_image_free(image_data);

        ImageMatrix original = original_uint.cast<double>();
        std::cout << "Task 1. Image size: " << original.rows() << " x " << original.cols() << " (rows x cols)\n";

        std::mt19937 rng(std::random_device{}());
        std::uniform_int_distribution<int> noise_dist(-40, 40);

        ImageMatrix noisy(original.rows(), original.cols());
        for (int r = 0; r < original.rows(); ++r) {
            for (int c = 0; c < original.cols(); ++c) {
                const double perturbed = original(r, c) + static_cast<double>(noise_dist(rng));
                noisy(r, c) = std::clamp(perturbed, 0.0, 255.0);
            }
        }

        ByteImageMatrix noisy_uint = noisy.unaryExpr([](double v) { return clamp_pixel(v); });
        save_png("Outputs/noise.png", noisy_uint);
        std::cout << "Task 2. Noise image saved to Outputs/noise.png\n";

        const Eigen::VectorXd v = Eigen::Map<const Eigen::VectorXd>(original.data(), original.size());
        const Eigen::VectorXd w = Eigen::Map<const Eigen::VectorXd>(noisy.data(), noisy.size());
        std::cout << std::setprecision(12);
        std::cout << "Task 3. v size = " << v.size() << ", ||v||_2 = " << v.norm() << "\n";

        const Eigen::MatrixXd hav1 = create_filter_hav1();
        const Eigen::SparseMatrix<double> A1 = create_convolution_matrix(hav1, original.rows(), original.cols());
        std::cout << "Task 4. Non-zero entries in A1 = " << A1.nonZeros() << "\n";

        Eigen::VectorXd smoothed_vec = A1 * w;
        ByteImageMatrix smoothed_img = vector_to_image(smoothed_vec, height, width);
        save_png("Outputs/smoothing.png", smoothed_img);
        std::cout << "Task 5. Smoothed image saved to Outputs/smoothing.png\n";

        const Eigen::MatrixXd hsh1 = create_filter_hsh1();
        const Eigen::SparseMatrix<double> A2 = create_convolution_matrix(hsh1, original.rows(), original.cols());
        const bool is_a2_symmetric = A2.isApprox(A2.transpose(), 1e-12);
        std::cout << "Task 6. Non-zero entries in A2 = " << A2.nonZeros()
                  << ", A2 symmetric? " << (is_a2_symmetric ? "true" : "false") << "\n";

        Eigen::VectorXd sharpened_vec = A2 * v;
        ByteImageMatrix sharpened_img = vector_to_image(sharpened_vec, height, width);
        save_png("Outputs/sharpening.png", sharpened_img);
        std::cout << "Task 7. Sharpened image saved to Outputs/sharpening.png\n";

        Eigen::saveMarket(A2, "Outputs/A2.mtx");
        export_vector_market("Outputs/w.mtx", w);
        std::cout << "Task 8. Exported A2.mtx and w.mtx under /Outputs\n";

        Eigen::VectorXd lis_solution;
        const LisSolveResult lis_result = solve_with_lis(A2, w, lis_solution);
        if (lis_result.success) {
            std::cout << "Task 8. LIS iterations = " << lis_result.iterations
                      << ", final residual = " << lis_result.residual << "\n";
        } else {
            std::cout << "Task 8. LIS solver failed (status code reported)." << '\n';
        }

        if (lis_result.success && lis_solution.size() == v.size()) {
            ByteImageMatrix lis_img = vector_to_image(lis_solution, height, width);
            save_png("Outputs/lis_solution.png", lis_img);
            std::cout << "Task 9. LIS solution image saved to Outputs/lis_solution.png\n";
        }

        const Eigen::MatrixXd hed2 = create_filter_hed2();
        const Eigen::SparseMatrix<double> A3 = create_convolution_matrix(hed2, original.rows(), original.cols());
        const bool is_a3_symmetric = A3.isApprox(A3.transpose(), 1e-12);
        std::cout << "Task 10. A3 symmetric? " << (is_a3_symmetric ? "true" : "false") << "\n";

        Eigen::VectorXd edges_vec = A3 * v;
        ByteImageMatrix edges_img = vector_to_image(edges_vec, height, width);
        save_png("Outputs/edge_detection.png", edges_img);
        std::cout << "Task 11. Edge detection image saved to Outputs/edge_detection.png\n";

        Eigen::SparseMatrix<double> system_matrix(A3.rows(), A3.cols());
        system_matrix.setIdentity();
        system_matrix *= 3.0;
        system_matrix += A3;

        Eigen::BiCGSTAB<Eigen::SparseMatrix<double>> eig_solver;
        eig_solver.setMaxIterations(1000);
        eig_solver.setTolerance(1e-8);
        eig_solver.compute(system_matrix);
        Eigen::VectorXd y = eig_solver.solve(w);
        if (eig_solver.info() == Eigen::Success) {
            std::cout << "Task 12. Eigen solver iterations = " << eig_solver.iterations()
                      << ", final residual = " << eig_solver.error() << "\n";
            ByteImageMatrix eigen_img = vector_to_image(y, height, width);
            save_png("Outputs/eigen_solution.png", eigen_img);
            std::cout << "Task 13. Eigen solution image saved to Outputs/eigen_solution.png\n";
        } else {
            std::cout << "Task 12. Eigen solver failed with status " << eig_solver.info() << "\n";
        }

        return 0;
    } catch (const std::exception &ex) {
        std::cerr << "Error: " << ex.what() << '\n';
        return 1;
    }
}
