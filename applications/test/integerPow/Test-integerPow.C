#include "cpuTime.H"
#include "Pair.H"
#include "Random.H"

using namespace Foam;

//- Number of calculations to do. Need 10^8-ish to get enough time (~2 seconds)
//  to meaningfully compare.
static const label n = 100000000;

//- Range of exponents to calculate. From -E to +E. Zero will be omitted so
//  that we can do roots, too. At small values (<8) integer powers win. At
//  larger values (>16), scalar powers do better.
static const label E = 4;

//- Run a given power method
template<class Pow>
Pair<scalar> run(Pow pow)
{
    Random rndGen(0);

    cpuTime time;

    scalar y = 0;

    for (label i = 0; i < n; ++ i)
    {
        const scalar x = rndGen.sampleAB<scalar>(rootSmall, rootGreat);
        const label e =
            rndGen.sample01<label>()
          ? rndGen.sampleAB<label>(- E, 0)
          : rndGen.sampleAB<label>(1, E + 1);

        y = y*scalar(i)/(i + 1) + pow(x, e)/(i + 1);
    }

    return Pair<scalar>(y, time.cpuTimeIncrement());
}

//- Scalar power function
scalar scalarPow(const scalar x, const scalar e)
{
    return Foam::pow(x, e);
}

//- Scalar root function
scalar scalarRoot(const scalar x, const scalar e)
{
    return Foam::pow(x, 1/e);
}

//- Run tests
int main()
{
    Info<< "Averaging " << n << " power/root calculations with exponents "
        << "ranging from " << - E << " to " << + E << ":" << endl;

    Info<< incrIndent;

    Info<< indent << " Scalar power: " << flush;
    const Pair<scalar> scalarPowYTime =
        run(static_cast<scalar (*)(scalar, scalar)>(&scalarPow));
    Info<< "Calculated value " << scalarPowYTime.first()
        << " in " << scalarPowYTime.second() << " s" << endl;

    Info<< indent << " Integer power: " << flush;
    const Pair<scalar> integerPowYTime =
        run(static_cast<scalar (*)(scalar, label)>(&integerPow));
    Info<< "Calculated value " << integerPowYTime.first()
        << " in " << integerPowYTime.second() << " s" << endl;

    Info<< indent << " Scalar root: " << flush;
    const Pair<scalar> scalarRootYTime =
        run(static_cast<scalar (*)(scalar, scalar)>(&scalarRoot));
    Info<< "Calculated value " << scalarRootYTime.first()
        << " in " << scalarRootYTime.second() << " s" << endl;

    Info<< indent << " Integer root: " << flush;
    const Pair<scalar> integerRootYTime =
        run(static_cast<scalar (*)(scalar, label)>(&integerRoot));
    Info<< "Calculated value " << integerRootYTime.first()
        << " in " << integerRootYTime.second() << " s" << endl;

    Info<< decrIndent;

    return 0;
}
