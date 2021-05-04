//
// Created by d-qql on 22.11.2020.
//

#ifndef PETGEN2020_MATRIXTRAITS_H
#define PETGEN2020_MATRIXTRAITS_H
#include <cstddef>
#include <type_traits>
#include "Product.h"
#include <vector>
#include <iostream>
#include <concepts>
#include <omp.h>
#include <ctime>

template<typename Matrix>
concept isApplyToVectorSpecified = requires (Matrix const &m){
    m.applyToVector({});
};

template<typename Matrix>
concept isAccessSpecified = requires (Matrix &m){
    m(0, 0);
};

template<typename Matrix>
concept isProduct = requires(Matrix const &m){
    m.left();
    m.right();
};
template<class Matrix>
struct MatrixTraits
{
    using elm_t = typename Matrix::elm_t; // Тип элементов матрицы
    using idx_t = typename Matrix::idx_t; // Тип индексов матрицы


    static auto sizeH(Matrix const &m)
    -> decltype(m.sizeH())
    {
        return m.sizeH();
    }

    static auto sizeW(Matrix const &m)
    -> decltype(m.sizeW())
    {
        return m.sizeW();
    }

    static elm_t& access(Matrix &m, const idx_t i,const  idx_t j) requires isAccessSpecified<Matrix>
    {
        return m(i, j);
    }

    [[maybe_unused]] static elm_t access(Matrix const &m, const idx_t i,const  idx_t j) requires isAccessSpecified<Matrix>
    {
        return m(i, j);
    }
    static bool isSqr(Matrix const &m){
        return m.sizeH()==m.sizeW();

    }
    static void add(Matrix &m, const idx_t i, const idx_t j, const elm_t value){
        m.add(i, j, value);
    }
    static void multiply(Matrix &m, const idx_t i, const idx_t j, const elm_t value){
        m.multiply(i, j, value);
    }

};


template<typename Matrix>
std::vector<typename MatrixTraits<Matrix>::elm_t> applyToVector(Matrix const &matrix, std::vector<typename MatrixTraits<Matrix>::elm_t> const &vec)
requires isProduct<Matrix> {
    return applyToVector(matrix.left(), applyToVector(matrix.right(), vec));
}
template<typename Matrix>
std::vector<typename MatrixTraits<Matrix>::elm_t> applyToVector(Matrix const &matrix, std::vector<typename MatrixTraits<Matrix>::elm_t> const &vec)
requires isApplyToVectorSpecified<Matrix>{
    return matrix.applyToVector(vec);
}
template<typename Matrix>
std::vector<typename MatrixTraits<Matrix>::elm_t> applyToVector(Matrix const &matrix, std::vector<typename MatrixTraits<Matrix>::elm_t> const &vec)
requires isAccessSpecified<Matrix> && (!isApplyToVectorSpecified<Matrix>){
    using elm_t = typename MatrixTraits<Matrix>::elm_t;
    using idx_t = typename MatrixTraits<Matrix>::idx_t;
    std::vector<elm_t> result;
    result.resize(matrix.sizeH());
    for(idx_t i = 0; i < matrix.sizeH(); ++i){
        auto c = static_cast<elm_t>(0);
        for(idx_t j = 0; j < matrix.sizeW(); ++j){
            c+=matrix(i, j) * vec[j];
        }
        result[i] = c;
    }
    return result;
}

template<typename Matrix>
[[maybe_unused]] auto trace(Matrix const &m)
{
    using MT = MatrixTraits<Matrix>;
    if(MT::isSqr(m)) {
        auto result = static_cast<typename MT::elm_t>(0);
        for (typename MT::idx_t i = 0; i < MT::sizeH(m); ++i)
            result += MT::access(m, i, i);
        return result;
    }else{
        std::cout<<"Can't calculate TRACE; Matrix is not square;"<<std::endl;
    }
}

template<typename Matrix>
[[maybe_unused]] void print(Matrix const &m){
    using MT = MatrixTraits<Matrix>;
    std::cout<<std::endl;
    for(typename MT::idx_t i = 0; i < MT::sizeH(m); ++i){
        for(typename MT::idx_t j = 0; j < MT::sizeW(m); ++j){
            std::cout<<MT::access(m, i, j)<<" ";
        }
        std::cout<<std::endl;
    }
}

template<typename Matrix, typename OtherMatrix>
Matrix& operator+=(Matrix& matrix, OtherMatrix const &other){
    using idx_t = typename MatrixTraits<Matrix>::idx_t;
    using MTcur = MatrixTraits<Matrix>;
    using MToth = MatrixTraits<OtherMatrix>;
    if(std::is_same<typename MTcur::elm_t, typename MToth::elm_t>::value && MToth::sizeH(other) == MTcur::sizeH(matrix) && MTcur::sizeW(matrix) == MToth::sizeW(other)){
        //using elm_t = typename MatrixTraits<Matrix>::elm_t;
        idx_t sizeH = MTcur::sizeH(matrix);
        idx_t sizeW = MTcur::sizeW(matrix);
        for(idx_t i = 0; i < sizeH; i++){
            for(idx_t j = 0; j < sizeW; j++){
                MTcur::add(matrix, i, j, MToth::access(other, i, j));
            }
        }
    }
    return matrix;
}

template<typename Matrix, typename OtherMatrix>
Matrix& operator-=(Matrix& matrix, OtherMatrix const &other){
    using idx_t = typename MatrixTraits<Matrix>::idx_t;
    using MTcur = MatrixTraits<Matrix>;
    using MToth = MatrixTraits<OtherMatrix>;
    if(std::is_same<typename MTcur::elm_t, typename MToth::elm_t>::value && MToth::sizeH(other) == MTcur::sizeH(matrix) && MTcur::sizeW(matrix) == MToth::sizeW(other)){
        //using elm_t = typename MatrixTraits<Matrix>::elm_t;
        idx_t sizeH = MTcur::sizeH(matrix);
        idx_t sizeW = MTcur::sizeW(matrix);
        for(idx_t i = 0; i < sizeH; i++){
            for(idx_t j = 0; j < sizeW; j++){
                MTcur::add(matrix, i, j, -MToth::access(other, i, j));
            }
        }
    }
    return matrix;
}
template<typename Matrix>
Matrix& operator*=(Matrix& matrix, typename MatrixTraits<Matrix>::elm_t k){
    using MT = MatrixTraits<Matrix>;
    using idx_t = typename MatrixTraits<Matrix>::idx_t;
    for (idx_t i = 0; i < MT::sizeH(); ++i) {
        for (idx_t j = 0; j < MT::sizeW(); ++j) {
            MT::multiply(matrix, i, j, k);
        }
    }
}

template<typename leftMatrix, typename rightMatrix>
Product<leftMatrix, rightMatrix> operator*(leftMatrix const &lft, rightMatrix const &rgh) {
    using MTL = MatrixTraits<leftMatrix>;
    using MTR = MatrixTraits<rightMatrix>;
    if (MTL::sizeW(lft) == MTR::sizeH(rgh)){
        Product<leftMatrix, rightMatrix> result = Product(lft, rgh);
        return result;
    }
}

template<typename T>
T operator*(std::vector<T> const &fvec, std::vector<T> const &svec){
    auto c = static_cast<T>(0);
    std::clock_t start, end;
    for(size_t i = 0; i < fvec.size(); ++i){
        c+=fvec[i]*svec[i];
    }
    return c;
}

template<typename T>
std::vector<T> operator-(std::vector<T> const &fvec, std::vector<T> const &svec){
    std::vector<T> result;
    result.resize(fvec.size());
    for(size_t i = 0; i<result.size(); ++i){
        result[i] = fvec[i] - svec[i];
    }
    return result;
}
template<typename T>
std::vector<T> operator+(std::vector<T> const &fvec, std::vector<T> const &svec){
    std::vector<T> result;
    result.resize(fvec.size());
    for(size_t i = 0; i<result.size(); ++i){
        result[i] = fvec[i] + svec[i];
    }
    return result;
}

template<typename T>
std::vector<T> operator*(T k, std::vector<T> const &vec){
    std::vector<T> result;
    result.resize(vec.size());
    for(size_t i = 0; i<result.size(); ++i){
        result[i] = k*vec[i];
    }
    return result;
}

template<typename T>
T norm(std::vector<T> const &x){
    T result = 0;
    for(int i = 0; i < x.size(); i++){
        result+=pow(x[i],2);
    }
    return sqrt(result);
}




#endif //PETGEN20_MATRIXTRAITS_H
