//
// Created by d-qql on 22.11.2020.
//

#ifndef PETGEN20_PRODUCT_H
#define PETGEN20_PRODUCT_H

#include <cmath>
#include <vector>
#include <concepts>



template<typename leftMatrix, typename rightMatrix>
class Product{
public:

    using idx_t = std::size_t;
    using elm_t = typename leftMatrix::elm_t;
    using LM = leftMatrix;
    using RM = rightMatrix;
private:
    leftMatrix L;
    rightMatrix R;
public:
    Product(LM l, RM r): L(l), R(r){

    }
    ~Product() = default;

    LM left() const{
        return L;
    }
    RM right() const{
        return R;
    }
    LM& left(){
        return L;
    }
    RM& right(){
        return R;
    }
    [[nodiscard]] idx_t sizeH() const{
        return L.sizeH();
    }

    [[nodiscard]] idx_t sizeW() const{
        return R.sizeW();
    }

};
#endif //PETGEN20_PRODUCT_H
