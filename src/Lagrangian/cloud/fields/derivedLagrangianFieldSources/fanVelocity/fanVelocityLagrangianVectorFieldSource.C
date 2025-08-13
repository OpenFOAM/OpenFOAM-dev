/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2025 OpenFOAM Foundation
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

#include "fanVelocityLagrangianVectorFieldSource.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fanVelocityLagrangianVectorFieldSource::
fanVelocityLagrangianVectorFieldSource
(
    const regIOobject& iIo,
    const dictionary& dict
)
:
    LagrangianVectorFieldSource(iIo, dict),
    Function1LagrangianFieldSource(*this),
    fanDirectionLagrangianVectorFieldSource(*this, dict),
    Ucentre_
    (
        Function1<vector>::New
        (
            "Ucentre",
            iIo.time().userUnits(),
            dimVelocity,
            dict
        )
    )
{}


Foam::fanVelocityLagrangianVectorFieldSource::
fanVelocityLagrangianVectorFieldSource
(
    const fanVelocityLagrangianVectorFieldSource& field,
    const regIOobject& iIo
)
:
    LagrangianVectorFieldSource(field, iIo),
    Function1LagrangianFieldSource(*this),
    fanDirectionLagrangianVectorFieldSource(field, *this),
    Ucentre_(field.Ucentre_, false)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fanVelocityLagrangianVectorFieldSource::
~fanVelocityLagrangianVectorFieldSource()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::LagrangianSubVectorField>
Foam::fanVelocityLagrangianVectorFieldSource::value
(
    const LagrangianInjection& injection,
    const LagrangianSubMesh& subMesh
) const
{
    const LagrangianSubVectorField Ucentre
    (
        value(injection, subMesh, Ucentre_())
    );

    const LagrangianSubScalarField magUcentre(mag(Ucentre));

    return magUcentre*direction(injection, Ucentre/magUcentre);
}


void Foam::fanVelocityLagrangianVectorFieldSource::write(Ostream& os) const
{
    LagrangianVectorFieldSource::write(os);

    fanDirectionLagrangianVectorFieldSource::write(os);

    writeEntry(os, db().time().userUnits(), dimVelocity, Ucentre_());
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makeLagrangianTypeFieldSource
    (
        LagrangianVectorFieldSource,
        fanVelocityLagrangianVectorFieldSource
    );
}

// ************************************************************************* //
