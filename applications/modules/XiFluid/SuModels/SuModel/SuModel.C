/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2024 OpenFOAM Foundation
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

#include "SuModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(SuModel, 0);
    defineRunTimeSelectionTable(SuModel, dictionary);
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

bool Foam::SuModel::readCoeffs(const dictionary&)
{
    return true;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::SuModel::SuModel
(
    const psiuMulticomponentThermo& thermo,
    const fluidThermoThermophysicalTransportModel& thermoTransport
)
:
    thermo_(thermo),
    thermoTransport_(thermoTransport),
    turbulence_(thermoTransport.momentumTransport()),
    Su_
    (
        IOobject
        (
            "Su",
            thermo_.mesh().time().name(),
            thermo_.mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        thermo_.mesh(),
        dimensionedScalar(dimVelocity, 0)
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::SuModel::~SuModel()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::SuModel::read(const dictionary& combustionProperties)
{
    return readCoeffs
    (
        combustionProperties.subDict("Su").optionalSubDict(type() + "Coeffs")
    );
}


// ************************************************************************* //
