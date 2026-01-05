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

#include "flowRateConeDiskVelocityLagrangianVectorFieldSource.H"
#include "diskInjection.H"
#include "flowRateNumberLagrangianScalarFieldSource.H"
#include "grouped.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::flowRateConeDiskVelocityLagrangianVectorFieldSource::
flowRateConeDiskVelocityLagrangianVectorFieldSource
(
    const regIOobject& iIo,
    const dictionary& dict
)
:
    LagrangianVectorFieldSource(iIo, dict),
    cloudLagrangianFieldSource(*this),
    coneDiskDirectionLagrangianVectorFieldSource(*this, dict)
{}


Foam::flowRateConeDiskVelocityLagrangianVectorFieldSource::
flowRateConeDiskVelocityLagrangianVectorFieldSource
(
    const flowRateConeDiskVelocityLagrangianVectorFieldSource& field,
    const regIOobject& iIo
)
:
    LagrangianVectorFieldSource(field, iIo),
    cloudLagrangianFieldSource(*this),
    coneDiskDirectionLagrangianVectorFieldSource(field, *this)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::flowRateConeDiskVelocityLagrangianVectorFieldSource::
~flowRateConeDiskVelocityLagrangianVectorFieldSource()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::LagrangianSubVectorField>
Foam::flowRateConeDiskVelocityLagrangianVectorFieldSource::value
(
    const LagrangianInjection& injection,
    const LagrangianSubMesh& subMesh
) const
{
    const flowRateNumberLagrangianScalarFieldSource& flowRateNumber =
        cloud<clouds::grouped>(injection)
       .number
       .sources()[injection.name()]
       .fieldSourceCast<flowRateNumberLagrangianScalarFieldSource>(injection);

    const Lagrangian::diskInjection& diskInjection =
        modelCast<Lagrangian::diskInjection>(injection);

    const LagrangianSubVectorField direction
    (
        this->direction(injection, subMesh)
    );

    return
        flowRateNumber.Q(injection, subMesh)
       /diskInjection.area()
       *direction
       /average(direction & diskInjection.axis());
}


void Foam::flowRateConeDiskVelocityLagrangianVectorFieldSource::write
(
    Ostream& os
) const
{
    LagrangianVectorFieldSource::write(os);

    coneDiskDirectionLagrangianVectorFieldSource::write(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makeLagrangianTypeFieldSource
    (
        LagrangianVectorFieldSource,
        flowRateConeDiskVelocityLagrangianVectorFieldSource
    );
}

// ************************************************************************* //
