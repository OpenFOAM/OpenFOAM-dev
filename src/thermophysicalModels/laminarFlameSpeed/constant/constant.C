/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2022 OpenFOAM Foundation
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

#include "constant.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace laminarFlameSpeedModels
{
    defineTypeNameAndDebug(constant, 0);

    addToRunTimeSelectionTable
    (
        laminarFlameSpeed,
        constant,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::laminarFlameSpeedModels::constant::constant
(
    const dictionary& dict,
    const psiuMulticomponentThermo& ct
)
:
    laminarFlameSpeed(dict, ct),

    Su_(dict.lookup("Su"))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::laminarFlameSpeedModels::constant::~constant()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::laminarFlameSpeedModels::constant::operator()() const
{
    return volScalarField::New
    (
        "Su0",
        psiuMulticomponentThermo_.T().mesh(),
        Su_
    );
}


// ************************************************************************* //
