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

#include "specificHeatCapacityLagrangianScalarFieldSource.H"
#include "basicLagrangianThermo.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::specificHeatCapacityLagrangianScalarFieldSource::
~specificHeatCapacityLagrangianScalarFieldSource()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::LagrangianSubScalarField>
Foam::specificHeatCapacityLagrangianScalarFieldSource::value
(
    const LagrangianInjection& injection,
    const LagrangianSubMesh& subMesh
) const
{
    const basicLagrangianThermo& thermo =
        db().lookupObject<basicLagrangianThermo>
        (
            IOobject::groupName
            (
                physicalProperties::typeName,
                internalGroup()
            )
        );

    return
        thermo.Cv
        (
            thermo.T().sources()[injection.name()].value(injection, subMesh),
            injection
        );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makeNullConstructableLagrangianTypeFieldSource
    (
        LagrangianScalarFieldSource,
        specificHeatCapacityLagrangianScalarFieldSource
    );
}

// ************************************************************************* //
