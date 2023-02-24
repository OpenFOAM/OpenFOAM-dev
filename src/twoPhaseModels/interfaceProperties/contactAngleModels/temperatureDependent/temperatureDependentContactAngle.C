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

#include "temperatureDependentContactAngle.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace contactAngleModels
{
    defineTypeNameAndDebug(temperatureDependent, 0);
    addToRunTimeSelectionTable
    (
        contactAngleModel,
        temperatureDependent,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::contactAngleModels::temperatureDependent::temperatureDependent
(
    const dictionary& dict
)
:
    TName_(dict.lookupOrDefault<word>("T", "T")),
    theta0_(Function1<scalar>::New("theta0", dict))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::contactAngleModels::temperatureDependent::~temperatureDependent()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::scalarField>
Foam::contactAngleModels::temperatureDependent::cosTheta
(
    const fvPatchVectorField& Up,
    const vectorField& nHat
) const
{
    return cos
    (
        degToRad
        (
            theta0_->value
            (
                Up.patch().lookupPatchField<volScalarField, scalar>(TName_)
            )
        )
    );
}


void Foam::contactAngleModels::temperatureDependent::write
(
    Ostream& os
) const
{
    writeEntry(os, theta0_());
}


// ************************************************************************* //
