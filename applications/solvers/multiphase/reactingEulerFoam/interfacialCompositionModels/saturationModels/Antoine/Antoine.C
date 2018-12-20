/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2015-2018 OpenFOAM Foundation
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

#include "Antoine.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace saturationModels
{
    defineTypeNameAndDebug(Antoine, 0);
    addToRunTimeSelectionTable(saturationModel, Antoine, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::saturationModels::Antoine::Antoine
(
    const dictionary& dict,
    const objectRegistry& db
)
:
    saturationModel(db),
    A_("A", dimless, dict),
    B_("B", dimTemperature, dict),
    C_("C", dimTemperature, dict)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::saturationModels::Antoine::~Antoine()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::saturationModels::Antoine::pSat
(
    const volScalarField& T
) const
{
    return
        dimensionedScalar(dimPressure, 1)
       *exp(A_ + B_/(C_ + T));
}


Foam::tmp<Foam::volScalarField>
Foam::saturationModels::Antoine::pSatPrime
(
    const volScalarField& T
) const
{
    return - pSat(T)*B_/sqr(C_ + T);
}


Foam::tmp<Foam::volScalarField>
Foam::saturationModels::Antoine::lnPSat
(
    const volScalarField& T
) const
{
    return A_ + B_/(C_ + T);
}


Foam::tmp<Foam::volScalarField>
Foam::saturationModels::Antoine::Tsat
(
    const volScalarField& p
) const
{
    return
        B_/(log(p*dimensionedScalar(dimless/dimPressure, 1)) - A_)
      - C_;
}


// ************************************************************************* //
