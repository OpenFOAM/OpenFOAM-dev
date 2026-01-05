/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2025-2026 OpenFOAM Foundation
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

#include "totalPressureVelocityLagrangianVectorFieldSource.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::totalPressureVelocityLagrangianVectorFieldSource::
totalPressureVelocityLagrangianVectorFieldSource
(
    const regIOobject& iIo,
    const dictionary& dict
)
:
    LagrangianVectorFieldSource(iIo, dict),
    cloudLagrangianFieldSource(*this),
    Function1LagrangianFieldSource(*this),
    totalPressureVelocityMagnitudeLagrangianScalarFieldSource
    (
        *this,
        *this,
        dict
    ),
    direction_
    (
        Function1<vector>::New
        (
            "direction",
            iIo.time().userUnits(),
            dimless,
            dict
        )
    )
{}


Foam::totalPressureVelocityLagrangianVectorFieldSource::
totalPressureVelocityLagrangianVectorFieldSource
(
    const totalPressureVelocityLagrangianVectorFieldSource& field,
    const regIOobject& iIo
)
:
    LagrangianVectorFieldSource(field, iIo),
    cloudLagrangianFieldSource(*this),
    Function1LagrangianFieldSource(*this),
    totalPressureVelocityMagnitudeLagrangianScalarFieldSource
    (
        field,
        *this,
        *this
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::totalPressureVelocityLagrangianVectorFieldSource::
~totalPressureVelocityLagrangianVectorFieldSource()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::LagrangianSubVectorField>
Foam::totalPressureVelocityLagrangianVectorFieldSource::value
(
    const LagrangianInjection& injection,
    const LagrangianSubMesh& subMesh
) const
{
    return
        Umag(injection, subMesh)
       *normalised(value(subMesh, dimless, direction_()));
}


void Foam::totalPressureVelocityLagrangianVectorFieldSource::write
(
    Ostream& os
) const
{
    LagrangianVectorFieldSource::write(os);

    totalPressureVelocityMagnitudeLagrangianScalarFieldSource::write(os);

    writeEntry(os, db().time().userUnits(), unitNone, direction_());
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makeLagrangianTypeFieldSource
    (
        LagrangianVectorFieldSource,
        totalPressureVelocityLagrangianVectorFieldSource
    );
}

// ************************************************************************* //
