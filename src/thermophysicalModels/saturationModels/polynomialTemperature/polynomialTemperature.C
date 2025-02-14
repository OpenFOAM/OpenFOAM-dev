/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2015-2025 OpenFOAM Foundation
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

#include "polynomialTemperature.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace saturationModels
{
    defineTypeNameAndDebug(polynomialTemperature, 0);
    addToRunTimeSelectionTable
    (
        saturationTemperatureModel,
        polynomialTemperature,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::saturationModels::polynomialTemperature::polynomialTemperature
(
    const dictionary& dict
)
:
    saturationTemperatureModel(),
    C_(dict.lookup("C<8>"))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::saturationModels::polynomialTemperature::~polynomialTemperature()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::scalarField>
Foam::saturationModels::polynomialTemperature::Tsat
(
    const scalarField& p
) const
{
    tmp<scalarField> tTsat(new scalarField(p.size(), scalar(0)));
    scalarField& Tsat = tTsat.ref();

    forAll(Tsat, celli)
    {
        Tsat[celli] = C_.value(p[celli]);
    }

    return tTsat;
}


Foam::tmp<Foam::scalarField>
Foam::saturationModels::polynomialTemperature::TsatPrime
(
    const scalarField& p
) const
{
    tmp<scalarField> tTsatPrime(new scalarField(p.size(), scalar(0)));
    scalarField& TsatPrime = tTsatPrime.ref();

    forAll(TsatPrime, celli)
    {
        TsatPrime[celli] = C_.derivative(p[celli]);
    }

    return tTsatPrime;
}


// ************************************************************************* //
