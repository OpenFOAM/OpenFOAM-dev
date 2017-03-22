#include <ctime>
#include <random>

#include "cubicEqn.H"
#include "IOstreams.H"
#include "stringList.H"

using namespace Foam;

scalar randomScalar(const scalar min, const scalar max)
{
    static_assert
    (
        sizeof(long) == sizeof(scalar),
        "Scalar and long are not the same size"
    );
    static std::default_random_engine generator(std::time(0));
    static std::uniform_int_distribution<long>
        distribution
        (
            std::numeric_limits<long>::min(),
            std::numeric_limits<long>::max()
        );
    scalar x;
    do
    {
        long i = distribution(generator);
        x = reinterpret_cast<scalar&>(i);
    }
    while (min > mag(x) || mag(x) > max || !std::isfinite(x));
    return x;
};

template <class Type>
void test(const Type& polynomialEqn, const scalar tol)
{
    Roots<Type::nComponents - 1> r = polynomialEqn.roots();

    const scalar nan = std::numeric_limits<scalar>::quiet_NaN();
    const scalar inf = std::numeric_limits<scalar>::infinity();

    FixedList<label, Type::nComponents - 1> t;
    FixedList<scalar, Type::nComponents - 1> v(nan);
    FixedList<scalar, Type::nComponents - 1> e(nan);
    bool ok = true;
    forAll(r, i)
    {
        t[i] = r.type(i);
        switch (t[i])
        {
            case roots::real:
                v[i] = polynomialEqn.value(r[i]);
                e[i] = polynomialEqn.error(r[i]);
                ok = ok && mag(v[i]) < tol*mag(e[i]);
                break;
            case roots::posInf:
                v[i] = + inf;
                e[i] = nan;
                break;
            case roots::negInf:
                v[i] = - inf;
                e[i] = nan;
                break;
            default:
                v[i] = e[i] = nan;
                break;
        }
    }

    if (!ok)
    {
        Info<< "Coeffs: " << polynomialEqn << endl
            << " Types: " << t << endl
            << " Roots: " << r << endl
            << "Values: " << v << endl
            << "Errors: " << e << endl << endl;
    }
}

int main()
{
    const int nTests = 1000000;
    for (int t = 0; t < nTests; ++ t)
    {
        test
        (
            cubicEqn
            (
                randomScalar(1e-50, 1e+50),
                randomScalar(1e-50, 1e+50),
                randomScalar(1e-50, 1e+50),
                randomScalar(1e-50, 1e+50)
            ),
            100
        );
    }

    Info << nTests << " cubics tested" << endl;

    return 0;
}
