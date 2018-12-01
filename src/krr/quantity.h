/*
 * quantity.h
 *
 *  Created on: 13.11.2016
 *      Author: klaus
 */

#ifndef QUANTITY_H_
#define QUANTITY_H_

#include <iostream>
#include <ratio>

namespace unit {

static constexpr int milli = 1e-3;

template<class T, class F, int D> class Quantity;

class Unit {

};

class Energy: public Unit {

};

class eV: public Energy {

};

class kcal: public Energy {

};

template<char ...C> struct UnitGenerator {
    static constexpr int D = 1;
    static Unit &create() {
        return NULL;
    }
};

template<> struct UnitGenerator<'k', 'c', 'a', 'l'> {
    typedef kcal type;
    static constexpr int D = 1;

    static Unit &create() {
        return *(new kcal());
    }
};

template<> struct UnitGenerator<'e', 'V'> {
    typedef eV type;
    static constexpr int D = 1;

    static Unit &create() {
        return *(new eV());
    }
};

template<class T, int D> struct QuantityGenerator {
    static Quantity<T, double, D> &create() {
        return *(new Quantity<T, double, D> );
    }
};

template<class T, class F = double, int D = 1> class Quantity: public std::ratio<
        1, D> {
    static_assert(std::is_base_of<Unit,T>::value, "Template type is not derived from class Unit");

    template<class I, int Pow> struct Pow10 {
        static constexpr int res = Pow10<I, Pow - 1>::res * 10;
    };

    template<class I> struct Pow10<I, 0> {
        static constexpr int res = 1;
    };

public:
    Quantity() :
            value(1.0) {
    }
    virtual ~Quantity() {
    }

    Quantity<T, F, D> &operator*(F scale) const {
        value *= scale;
        return *this;
    }

    template<class F2, int D2> Quantity<T, F, D> &operator*(
            const Quantity<T, F2, D2> &q) const {
        Quantity<T, double, D * D2> *nq = new Quantity<T, double, D * D2>;
        nq.value = value * Pow10<double, D - D2>::res;
        return *nq;
    }

    template<class F2, int D2> friend Quantity<T, F, D> operator*(
            const double scalar, const Quantity<T, F2, D2>& quan) {
        Quantity<T, F2, D2> res;
        res.value *= scalar;
        return res;
    }

    const F operator()() const {
        return value;
    }

private:
    F value;
};

}

template<class T, T ... Chars>
constexpr unit::Quantity<typename unit::UnitGenerator<Chars...>::type> &operator""_u() {
    return unit::QuantityGenerator<typename unit::UnitGenerator<Chars...>::type,
            unit::UnitGenerator<Chars...>::D>::create();
}

#endif /* QUANTITY_H_ */
