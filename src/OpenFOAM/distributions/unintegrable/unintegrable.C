/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2023-2024 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "unintegrable.H"
#include "quadraticEqn.H"
#include "SubField.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

Foam::tmp<Foam::scalarField> Foam::distributions::unintegrable::integrate
(
    const scalarField& x,
    const scalarField& y
)
{
    tmp<scalarField> tResult(new scalarField(x.size()));
    scalarField& result = tResult.ref();

    result[0] = 0;

    for (label i = 1; i < x.size(); ++ i)
    {
        result[i] = result[i - 1] + (x[i] - x[i - 1])*(y[i] + y[i - 1])/2;
    }

    return tResult;
}


Foam::tmp<Foam::scalarField> Foam::distributions::unintegrable::integrateX
(
    const scalarField& x,
    const scalarField& y
)
{
    tmp<scalarField> tResult(new scalarField(x.size()));
    scalarField& result = tResult.ref();

    result[0] = 0;

    for (label i = 1; i < x.size(); ++ i)
    {
        const scalar x0 = x[i - 1], x1 = x[i];
        const scalar y0 = y[i - 1], y1 = y[i];

        result[i] =
            result[i - 1] + (x1 - x0)*(2*x0*y0 + x0*y1 + x1*y0 + 2*x1*y1)/6;
    }

    return tResult;
}


Foam::tmp<Foam::scalarField>
Foam::distributions::unintegrable::interpolateIntegrateXPow
(
    const scalarField& xStar,
    const label e,
    const scalarField& yStar,
    const scalarField& x
)
{
    tmp<scalarField> tResult(new scalarField(x.size()));
    scalarField& result = tResult.ref();

    label i = 0;

    while (i < x.size() && x[i] < xStar[0])
    {
        result[i] = 0;
        i ++;
    }

    scalar integral_PDFxPowE_0_j = 0;

    for (label iStar = 0; iStar < xStar.size() - 1; ++ iStar)
    {
        const scalar xPowE1_j = integerPow(xStar[iStar], e + 1);

        auto integral_xPowE1_j_x = [&](const scalar x)
        {
            const scalar xPowE1_i = integerPow(x, e + 1);

            const scalar integral_xPowE_j_x =
                e + 1 == 0
              ? log(x/xStar[iStar])
              : (xPowE1_i - xPowE1_j)/(e + 1);
            const scalar integral_xPowE1_j_x =
                e + 2 == 0
              ? log(x/xStar[iStar])
              : (xPowE1_i*x - xPowE1_j*xStar[iStar])/(e + 2);

            return
                yStar[iStar]*integral_xPowE_j_x
              + (yStar[iStar + 1] - yStar[iStar])
               /(xStar[iStar + 1] - xStar[iStar])
               *(integral_xPowE1_j_x - xStar[iStar]*integral_xPowE_j_x);
        };

        while (i < x.size() && x[i] < xStar[iStar + 1])
        {
            result[i] = integral_PDFxPowE_0_j + integral_xPowE1_j_x(x[i]);

            i ++;
        }

        integral_PDFxPowE_0_j += integral_xPowE1_j_x(xStar[iStar + 1]);
    }

    while (i < x.size())
    {
        result[i] = integral_PDFxPowE_0_j;
        i ++;
    }

    return tResult;
}


Foam::scalar Foam::distributions::unintegrable::sampleInterval
(
    const Pair<scalar>& x,
    const Pair<scalar>& Phi,
    const scalar s
)
{
    // Do a linear interpolation
    scalar f = (s - Phi[0])/(Phi[1] - Phi[0]);

    // Interpolate
    return (1 - f)*x[0] + f*x[1];
}


Foam::scalar Foam::distributions::unintegrable::sampleInterval
(
    const Pair<scalar>& x,
    const Pair<scalar>& phi,
    const Pair<scalar>& Phi,
    const scalar s
)
{
    // Do a linear interpolation
    scalar f = (s - Phi[0])/(Phi[1] - Phi[0]);

    // Attempt to improve by doing a quadratic interpolation
    const scalar a = (x[1] - x[0])*(phi[1] - phi[0])/2;
    const scalar b = Phi[1] - Phi[0] - a;
    const scalar c = Phi[0] - s;
    const quadraticEqn eqn(a, b, c);
    const Roots<2> roots = eqn.roots();
    forAll(roots, rooti)
    {
        if
        (
            roots.type(rooti) == rootType::real
         && roots[rooti] > 0
         && roots[rooti] < 1
        )
        {
            f = roots[rooti];
        }
    }

    // Interpolate
    return (1 - f)*x[0] + f*x[1];
}


Foam::scalar Foam::distributions::unintegrable::sample
(
    const scalarField& x,
    const scalarField& Phi,
    const scalar s
)
{
    // Find the interval
    label i = 0;
    for (; i < x.size() && Phi[i + 1] < s; ++ i);

    // Sample
    return
        sampleInterval
        (
            {x[i], x[i + 1]},
            {Phi[i], Phi[i + 1]},
            s
        );
}


Foam::scalar Foam::distributions::unintegrable::sample
(
    const scalarField& x,
    const scalarField& phi,
    const scalarField& Phi,
    const scalar s
)
{
    // Find the interval
    label i = 0;
    for (; i < x.size() && Phi[i + 1] < s; ++ i);

    // Sample
    return
        sampleInterval
        (
            {x[i], x[i + 1]},
            {phi[i], phi[i + 1]},
            {Phi[i], Phi[i + 1]},
            s
        );
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

const Foam::scalarField& Foam::distributions::unintegrable::x() const
{
    if (xPtr_.valid()) return xPtr_();

    xPtr_.set(new scalarField(n_));
    scalarField& x = xPtr_();

    forAll(x, i)
    {
        const scalar f = scalar(i)/(n_ + 1);
        x[i] = (1 - f)*this->min() + f*this->max();
    }

    static const scalar tol = sqrt(sqrt(pow3(small)));

    for (label iter = 0; iter < ceil(std::log2(1/tol)); ++ iter)
    {
        const scalarField phi(this->phi(this->q(), x));
        const scalarField Phi(integrate(x, phi));

        const scalar dPhi = (Phi.last() - Phi.first())/(n_ - 1);

        if (distribution::debug)
        {
            scalar error = -vGreat;
            for (label i = 0; i < n_ - 1; ++ i)
            {
                error = Foam::max(error, (1 - (Phi[i + 1] - Phi[i])/dPhi)/n_);
            }

            Info<< indent << "Interval spacing iteration #" << iter
                << ", error=" << error << endl;
        }

        const scalarField xPrev(x);

        x.first() = this->min();

        for (label i0 = 0, i = 1; i0 < n_ - 1; ++ i0)
        {
            while (Phi[i0 + 1] > i*dPhi)
            {
                const scalar xNext =
                    sampleInterval
                    (
                        {xPrev[i0], xPrev[i0 + 1]},
                        {phi[i0], phi[i0 + 1]},
                        {Phi[i0], Phi[i0 + 1]},
                        i*dPhi
                    );

                x[i] = (x[i] + xNext)/2;

                i ++;
            }
        }

        x.last() = this->max();
    }

    return x;
}


const Foam::scalarField& Foam::distributions::unintegrable::PDF() const
{
    if (PDFPtr_.valid()) return PDFPtr_();

    const Pair<scalar> Phi01 = this->Phi01();

    PDFPtr_.set((phi(this->q(), x())/(Phi01[1] - Phi01[0])).ptr());

    return PDFPtr_();
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

Foam::tmp<Foam::scalarField> Foam::distributions::unintegrable::Phi
(
    const label q,
    const scalarField& x
) const
{
    return integrate(x, phi(q, x));
}


Foam::Pair<Foam::scalar> Foam::distributions::unintegrable::Phi01
(
    const label q
) const
{
    const scalarField Phi(this->Phi(this->q(), x()));
    return Pair<scalar>(Phi.first(), Phi.last());
}


const Foam::Pair<Foam::scalar>&
Foam::distributions::unintegrable::Phi01() const
{
    if (Phi01Ptr_.valid()) return Phi01Ptr_();

    Phi01Ptr_.set(new Pair<scalar>(Phi01(this->q())));

    return Phi01Ptr_();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::distributions::unintegrable::unintegrable
(
    const word& name,
    const unitConversion& units,
    const dictionary& dict,
    const label sampleQ,
    randomGenerator&& rndGen
)
:
    distribution(name, units, dict, sampleQ, std::move(rndGen)),
    n_((1<<dict.lookupOrDefault<label>("level", 16)) + 1)
{}


Foam::distributions::unintegrable::unintegrable
(
    const label Q,
    const label sampleQ,
    randomGenerator&& rndGen,
    const label n
)
:
    distribution(Q, sampleQ, std::move(rndGen)),
    n_(n)
{}


Foam::distributions::unintegrable::unintegrable
(
    const unintegrable& d,
    const label sampleQ
)
:
    distribution(d, sampleQ),
    n_(d.n_)
{
    // Copy the data over if it exists and if the sampling moment is the same,
    // otherwise leave it to be re-generated
    if (q() == d.q())
    {
        if (d.xPtr_.valid())
        {
            xPtr_.set(new scalarField(d.xPtr_()));
        }
        if (d.Phi01Ptr_.valid())
        {
            Phi01Ptr_.set(new Pair<scalar>(d.Phi01Ptr_()));
        }
        if (d.PDFPtr_.valid())
        {
            PDFPtr_.set(new scalarField(d.PDFPtr_()));
        }
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::distributions::unintegrable::~unintegrable()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::distributions::unintegrable::sample() const
{
    const scalar s = this->rndGen_.template sample01<scalar>();

    const scalar dCDF = scalar(1)/(n_ - 1);
    const label samplei = floor(s/dCDF);

    return
        sampleInterval
        (
            {x()[samplei], x()[samplei + 1]},
            {PDF()[samplei], PDF()[samplei + 1]},
            {samplei*dCDF, (samplei + 1)*dCDF},
            s
        );
}


Foam::scalar Foam::distributions::unintegrable::mean() const
{
    const Pair<scalar> Phi01 = this->Phi01();
    const scalarField Mu(Phi(this->q() + 1, x()));
    const Pair<scalar> Mu01(Mu.first(), Mu.last());
    return (Mu01[1] - Mu01[0])/(Phi01[1] - Phi01[0]);
}


Foam::tmp<Foam::scalarField>
Foam::distributions::unintegrable::integralPDFxPow
(
    const scalarField& x,
    const label e,
    const bool
) const
{
    return interpolateIntegrateXPow(this->x(), e, this->PDF(), x);
}


void Foam::distributions::unintegrable::write
(
    Ostream& os,
    const unitConversion& units
) const
{
    distribution::write(os, units);

    // Recover the level
    label n = n_ - 1, level = 0;
    while (n >>= 1) ++ level;

    writeEntryIfDifferent(os, "level", label(16), level);
}


Foam::tmp<Foam::scalarField>
Foam::distributions::unintegrable::plotPDF(const scalarField& x) const
{
    const scalarField phi(this->phi(this->q(), x));
    const Pair<scalar> Phi01 = this->Phi01();
    return this->clipPDF(x, phi/(Phi01[1] - Phi01[0]));
}


// ************************************************************************* //
