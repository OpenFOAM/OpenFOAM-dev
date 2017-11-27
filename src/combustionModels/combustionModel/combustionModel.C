/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2017 OpenFOAM Foundation
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
}

const Foam::word Foam::combustionModel::combustionPropertiesName
(
    "combustionProperties"
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::combustionModel::combustionModel
(
    const word& modelType,
    basicThermo& thermo,
    const compressibleTurbulenceModel& turb,
    const word& combustionProperties
)
:
    IOdictionary
    (
        IOobject
        (
            thermo.phasePropertyName(combustionProperties),
            thermo.db().time().constant(),
            thermo.db(),
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    ),
    mesh_(thermo.p().mesh()),
    turb_(turb),
    active_(lookupOrDefault<Switch>("active", true)),
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
        this->lookup("active") >> active_;
        coeffs_ = optionalSubDict(modelType_ + "Coeffs");
        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
