/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2026 OpenFOAM Foundation
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

#include "pressureLagrangianScalarFieldSource.H"
#include "fluidLagrangianThermo.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::pressureLagrangianScalarFieldSource::
~pressureLagrangianScalarFieldSource()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::LagrangianSubScalarField>
Foam::pressureLagrangianScalarFieldSource::value
(
    const LagrangianInjection& injection,
    const LagrangianSubMesh& subMesh
) const
{
    return
        db().lookupObject<fluidLagrangianThermo>
        (
            IOobject::groupName
            (
                physicalProperties::typeName,
                internalGroup()
            )
        ).p(injection, subMesh);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makeNullConstructableLagrangianTypeFieldSource
    (
        LagrangianScalarFieldSource,
        pressureLagrangianScalarFieldSource
    );
}

// ************************************************************************* //
