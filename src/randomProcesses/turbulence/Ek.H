#ifndef Ek_H
#define Ek_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

inline tmp<scalarField> Ek
(
    const scalar Ea,
    const scalar k0,
    const scalarField& k
)
{
    tmp<scalarField> tEk = Ea*pow(k/k0, 4.0)*exp(-2.0*sqr(k/k0));

    /*
    scalarField& Ekf = tEk();

    label i;
    forAll(Ekf, i)
    {
        if (k[i] < 2 || k[i] > 10)
        {
            Ekf[i] = 0.0;
        }
    }
    */

    return tEk;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#endif
