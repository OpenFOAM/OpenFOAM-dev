#include <iostream>

template<class C>
class Vector;

template<class C>
Vector<C> operator+(const Vector<C>& v1, const Vector<C>& v2);

template<class C>
std::ostream& operator<<(std::ostream& os, const Vector<C>& v);


/*---------------------------------------------------------------------------*\
                           Class Vector Declaration
\*---------------------------------------------------------------------------*/

template<class C>
class Vector
{

    double X, Y;

public:

    inline Vector(const double x, const double y);

    C x() const
    {
        return X;
    }

    C y() const
    {
        return Y;
    }

    friend Vector<C> operator+ <C>(const Vector<C>& v1, const Vector<C>& v2);

    friend std::ostream& operator<<(std::ostream& os, const Vector<C>& v)
    {
        os  << v.X << '\t' << v.Y << '\n';
        return os;
    }
};

template<class C>
inline Vector<C>::Vector(const double x, const double y)
{
    X = x;
    Y = y;
}


template<class C>
inline Vector<C> operator+(const Vector<C>& v1, const Vector<C>& v2)
{
    return Vector<C>(v1.X+v2.X, v1.Y+v2.Y);
}
