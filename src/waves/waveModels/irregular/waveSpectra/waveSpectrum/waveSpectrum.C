/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2023 OpenFOAM Foundation
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

#include "waveSpectrum.H"
#include "unintegrable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(waveSpectrum, 0);
    defineRunTimeSelectionTable(waveSpectrum, dictionary);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::tmp<Foam::scalarField> Foam::waveSpectrum::refined
(
    const label n,
    const scalarField& x
)
{
    if (n == x.size())
    {
        return x;
    }

    tmp<scalarField> txx(new scalarField(n, -vGreat));
    scalarField& xx = txx.ref();

    if (n > x.size())
    {
        for (label i = 0; i < x.size() - 1; ++ i)
        {
            const label j0 = i*(n - 1)/(x.size() - 1);
            const label j1 = (i + 1)*(n - 1)/(x.size() - 1);

            xx[j0] = x[i];

            for (label j = j0 + 1; j < j1; ++ j)
            {
                const scalar f = scalar(j - j0)/(j1 - j0);

                xx[j] = (1 - f)*x[i] + f*x[i + 1];
            }
        }

        xx.last() = x.last();
    }
    else
    {
        forAll(xx, j)
        {
            const label i = j*(x.size() - 1)/(n - 1);

            xx[j] = x[i];
        }
    }

    return txx;
}


Foam::tmp<Foam::scalarField> Foam::waveSpectrum::refined
(
    const label n,
    const tmp<scalarField>& x
)
{
    if (n == x().size())
    {
        return x;
    }

    return refined(n, x());
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

Foam::scalar Foam::waveSpectrum::fFraction
(
    const scalar fraction,
    const scalar f1
) const
{
    // Integrate from zero to the maximum frequency
    const scalarField ff(refined(n_, scalarField(scalarList({small*f1, f1}))));
    const scalarField integralSS(integralS(ff));

    // Sample the integrated values for the frequency at the given fraction
    return
        distributions::unintegrable::sample
        (
            ff,
            integralSS/integralSS.last(),
            fraction
        );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::waveSpectrum::waveSpectrum(const waveSpectrum& spectrum)
:
    g_(spectrum.g_),
    n_(spectrum.n_)
{}


Foam::waveSpectrum::waveSpectrum(const dictionary& dict, const scalar g)
:
    g_(g),
    n_((1<<dict.lookupOrDefault<label>("level", 16)) + 1)
{}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::waveSpectrum>
Foam::waveSpectrum::New(const dictionary& dict, const scalar g)
{
    const word type = dict.lookup<word>("spectrum");

    if (debug)
    {
        Info<< "Selecting " << waveSpectrum::typeName << " " << type << endl;
    }

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(type);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown " << waveSpectrum::typeName << " " << type << nl << nl
            << "Valid spectrum types are:" << nl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return cstrIter()(dict.optionalSubDict(type + "Coeffs"), g);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::waveSpectrum::~waveSpectrum()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::scalarField> Foam::waveSpectrum::integralS
(
    const scalarField& f
) const
{
    const scalarField ff(refined(n_, f));
    return
        refined
        (
            f.size(),
            distributions::unintegrable::integrate(ff, S(ff))
        );
}


Foam::tmp<Foam::scalarField> Foam::waveSpectrum::integralFS
(
    const scalarField& f
) const
{
    const scalarField ff(refined(n_, f));
    return
        refined
        (
            f.size(),
            distributions::unintegrable::integrateX(ff, S(ff))
        );
}


void Foam::waveSpectrum::write(Ostream& os) const
{}


// ************************************************************************* //
