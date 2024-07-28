#ifndef AUT_AP_2024_Spring_HW1
#define AUT_AP_2024_Spring_HW1

// Matrix data structure
#include <vector>
#include <optional>
#include <random>
#include <string>
#include <stdexcept>

namespace algebra {
    template<typename T>
    using MATRIX = std::vector<std::vector<T>>;

    // Matrix initialization types
    enum class MatrixType { Zeros, Ones, Identity, Random };

    // 對整數型別的特化
    template<typename T, typename std::enable_if<std::is_integral<T>::value, int>::type = 0>
    T randm(T low, T high) {
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_int_distribution<T> dist(low, high);
        return dist(gen);
    }

    // 對浮點型別的特化
    template<typename T, typename std::enable_if<std::is_floating_point<T>::value, int>::type = 0>
    T randm(T low, T high) {
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<T> dist(low, high);
        return dist(gen);
    }

    // 主模板，處理 std::optional 的情況
    template<typename T>
    T randm(std::optional<T>& low, std::optional<T>& high) {
        if (low && high) {
            if constexpr (std::is_integral<T>::value) {
                std::random_device rd;
                std::mt19937 gen(rd());
                std::uniform_int_distribution<T> dist(*low, *high);
                return dist(gen);
            } else if constexpr (std::is_floating_point<T>::value) {
                std::random_device rd;
                std::mt19937 gen(rd());
                std::uniform_real_distribution<T> dist(*low, *high);
                return dist(gen);
            }
        }
        throw std::invalid_argument("Bounds are not set properly.");
    }

    // Function template for matrix initialization
    template<typename T>
    MATRIX<T> create_matrix(std::size_t rows, std::size_t columns, std::optional<MatrixType> type = MatrixType::Zeros,
                            std::optional<T> lowerBound = std::nullopt, std::optional<T> upperBound = std::nullopt);
    // sum_sub
    template<typename T>
    MATRIX<T> sum_sub(const MATRIX<T>& matrixA, const MATRIX<T>& matrixB, std::optional<std::string> operation = "sum");

    template<typename T>
    MATRIX<T> multiply(const MATRIX<T>& matrix, const T scalar);

    template<typename T>
    MATRIX<T> multiply(const MATRIX<T>& matrixA, const MATRIX<T>& matrixB);

    template<typename T>
    MATRIX<T> hadamard_product(const MATRIX<T>& matrixA, const MATRIX<T>& matrixB);

    template<typename T>
    MATRIX<T> transpose(const MATRIX<T>& matrix);

    template<typename T>
    T trace(const MATRIX<T>& matrix);

    template<typename T>
    double determinant(const MATRIX<T>& matrix);

    template<typename T>
    MATRIX<double> inverse(const MATRIX<T>& matrix);

    // display
    template<typename T>
    void display(const MATRIX<T>& matrix);
}
#endif //AUT_AP_2024_Spring_HW1
