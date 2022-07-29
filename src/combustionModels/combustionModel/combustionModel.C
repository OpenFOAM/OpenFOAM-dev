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

#include "combustionModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(combustionModel, 0);
    defineRunTimeSelectionTable(combustionModel, dictionary)
}

const Foam::word Foam::combustionModel::combustionPropertiesName
(
    "combustionProperties"
);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::IOobject Foam::combustionModel::createIOobject
(
    const fluidMulticomponentThermo& thermo,
    const word& combustionProperties
) const
{
    typeIOobject<IOdictionary> io
    (
        thermo.phasePropertyName(combustionProperties),
        thermo.T().mesh().time().constant(),
        thermo.T().mesh(),
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    );

    if (io.headerOk())
    {
        io.readOpt() = IOobject::MUST_READ_IF_MODIFIED;
        return io;
    }
    else
    {
        io.readOpt() = IOobject::NO_READ;
        return io;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::combustionModel::combustionModel
(
    const word& modelType,
    const fluidMulticomponentThermo& thermo,
    const compressibleMomentumTransportModel& turb,
    const word& combustionProperties
)
:
    IOdictionary(createIOobject(thermo, combustionProperties)),
    mesh_(thermo.T().mesh()),
    thermo_(thermo),
    turb_(turb),
    coeffs_(optionalSubDict(modelType + "Coeffs")),
    modelType_(modelType)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::combustionModel::~combustionModel()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::combustionModel::read()
{
    if (regIOobject::read())
    {
        coeffs_ = optionalSubDict(modelType_ + "Coeffs");
        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
