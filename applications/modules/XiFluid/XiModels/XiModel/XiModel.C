/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2024 OpenFOAM Foundation
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

#include "XiModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(XiModel, 0);
    defineRunTimeSelectionTable(XiModel, dictionary);
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

bool Foam::XiModel::readCoeffs(const dictionary&)
{
    return true;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::XiModel::XiModel
(
    const psiuMulticomponentThermo& thermo,
    const fluidThermoThermophysicalTransportModel& thermoTransport,
    const volScalarField& Su
)
:
    thermo_(thermo),
    thermoTransport_(thermoTransport),
    turbulence_(thermoTransport.momentumTransport()),
    Su_(Su),
    rho_(turbulence_.rho()),
    b_(thermo_.Y("b")),
    Xi_
    (
        IOobject
        (
            "Xi",
            thermo_.mesh().time().name(),
            thermo_.mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        thermo_.mesh(),
        dimensionedScalar("1", dimless, 1)
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::XiModel::~XiModel()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::XiModel::read(const dictionary& combustionProperties)
{
    return readCoeffs
    (
        combustionProperties.subDict("Xi").optionalSubDict(type() + "Coeffs")
    );
}


void Foam::XiModel::reset()
{
    for (label i=0; i<=Xi_.nOldTimes(); i++)
    {
        Xi_.oldTimeRef(i) = 1;
    }
}


// ************************************************************************* //
