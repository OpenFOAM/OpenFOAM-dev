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

#include "coneDiskVelocityLagrangianVectorFieldSource.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::coneDiskVelocityLagrangianVectorFieldSource::
coneDiskVelocityLagrangianVectorFieldSource
(
    const regIOobject& iIo,
    const dictionary& dict
)
:
    LagrangianVectorFieldSource(iIo, dict),
    Function1LagrangianFieldSource(*this),
    coneDiskDirectionLagrangianVectorFieldSource(*this, dict),
    Umag_
    (
        Function1<scalar>::New
        (
            "Umag",
            iIo.time().userUnits(),
            dimVelocity,
            dict
        )
    )
{}


Foam::coneDiskVelocityLagrangianVectorFieldSource::
coneDiskVelocityLagrangianVectorFieldSource
(
    const coneDiskVelocityLagrangianVectorFieldSource& field,
    const regIOobject& iIo
)
:
    LagrangianVectorFieldSource(field, iIo),
    Function1LagrangianFieldSource(*this),
    coneDiskDirectionLagrangianVectorFieldSource(field, *this),
    Umag_(field.Umag_, false)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::coneDiskVelocityLagrangianVectorFieldSource::
~coneDiskVelocityLagrangianVectorFieldSource()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::LagrangianSubVectorField>
Foam::coneDiskVelocityLagrangianVectorFieldSource::value
(
    const LagrangianInjection& injection,
    const LagrangianSubMesh& subMesh
) const
{
    return
        value(injection, subMesh, dimVelocity, Umag_())
       *direction(injection, subMesh);
}


void Foam::coneDiskVelocityLagrangianVectorFieldSource::write(Ostream& os) const
{
    LagrangianVectorFieldSource::write(os);

    coneDiskDirectionLagrangianVectorFieldSource::write(os);

    writeEntry(os, db().time().userUnits(), dimVelocity, Umag_());
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makeLagrangianTypeFieldSource
    (
        LagrangianVectorFieldSource,
        coneDiskVelocityLagrangianVectorFieldSource
    );
}

// ************************************************************************* //
