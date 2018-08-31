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
    static std::uniform_int_distribution<long> distribution
    (
        std::numeric_limits<long>::min(),
        std::numeric_limits<long>::max()
    );
    scalar x;
    do
    {
        long i = distribution(generator);
        scalar* ptr = reinterpret_cast<scalar*>(&i);
        x = *ptr;
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

    FixedList<rootType, Type::nComponents - 1> t;
    FixedList<scalar, Type::nComponents - 1> v(nan);
    FixedList<scalar, Type::nComponents - 1> e(nan);
    bool ok = true;
    forAll(r, i)
    {
        t[i] = r.type(i);
        switch (t[i])
        {
            case rootType::real:
                v[i] = polynomialEqn.value(r[i]);
                e[i] = polynomialEqn.error(r[i]);
                ok = ok && mag(v[i]) <= tol*mag(e[i]);
                break;
            case rootType::posInf:
                v[i] = + inf;
                e[i] = nan;
                break;
            case rootType::negInf:
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
         // << " Types: " << t << endl
            << " Roots: " << r << endl
            << "Values: " << v << endl
            << "Errors: " << e << endl << endl;
    }
}

int main()
{
    const scalar tol = 5;

    const label nTests = 1000000;
    for (label t = 0; t < nTests; ++ t)
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
            tol
        );
    }
    Info << nTests << " random cubics tested" << endl;

    const label coeffMin = -9, coeffMax = 10, nCoeff = coeffMax - coeffMin;
    for (label a = coeffMin; a < coeffMax; ++ a)
    {
        for (label b = coeffMin; b < coeffMax; ++ b)
        {
            for (label c = coeffMin; c < coeffMax; ++ c)
            {
                for (label d = coeffMin; d < coeffMax; ++ d)
                {
                    test(cubicEqn(a, b, c, d), tol);
                }
            }
        }
    }
    Info<< nCoeff*nCoeff*nCoeff*nCoeff << " integer cubics tested" << endl;

    return 0;
}
