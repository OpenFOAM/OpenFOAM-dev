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
#include "diskInjection.H"
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
    Function1LagrangianFieldSource<vector>(*this),
    Umag_
    (
        Function1<scalar>::New
        (
            "Umag",
            iIo.time().userUnits(),
            dimVelocity,
            dict
        )
    ),
    thetaInner_
    (
        Function1<scalar>::New
        (
            "thetaInner",
            iIo.time().userUnits(),
            unitDegrees,
            dict
        )
    ),
    thetaOuter_
    (
        Function1<scalar>::New
        (
            "thetaOuter",
            iIo.time().userUnits(),
            unitDegrees,
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
    Function1LagrangianFieldSource<vector>(*this),
    Umag_(field.Umag_, false),
    thetaInner_(field.thetaInner_, false),
    thetaOuter_(field.thetaOuter_, false)
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
    const Lagrangian::diskInjection& diskInjection =
        modelCast<Lagrangian::diskInjection>(injection);

    // Evaluate the axial velocity
    const LagrangianSubScalarField Umag
    (
        value(injection, subMesh, dimVelocity, Umag_())
    );

    // Get the geometry from the disk injection model
    const LagrangianSubScalarField rFrac
    (
        LagrangianSubScalarField::New
        (
            "r",
            subMesh,
            dimless,
            diskInjection.rFrac()
        )
    );
    const LagrangianSubVectorField axis
    (
        LagrangianSubVectorField::New
        (
            "axis",
            subMesh,
            dimless,
            diskInjection.axis()
        )
    );
    const LagrangianSubVectorField radial
    (
        LagrangianSubVectorField::New
        (
            "radial",
            subMesh,
            dimless,
            diskInjection.radial()
        )
    );

    // Evaluate the cone angle
    const tmp<LagrangianSubScalarField> tthetaInner =
        value(injection, subMesh, dimless, thetaInner_());
    const tmp<LagrangianSubScalarField> tthetaOuter =
        value(injection, subMesh, dimless, thetaOuter_());
    const LagrangianSubScalarField theta
    (
        (1 - rFrac)*tthetaInner + rFrac*tthetaOuter
    );

    // Return the velocity in the calculated direction
    return Umag*(cos(theta)*axis + sin(theta)*radial);
}


void Foam::coneDiskVelocityLagrangianVectorFieldSource::write(Ostream& os) const
{
    LagrangianVectorFieldSource::write(os);

    writeEntry(os, db().time().userUnits(), dimVelocity, Umag_());
    writeEntry(os, db().time().userUnits(), unitDegrees, thetaInner_());
    writeEntry(os, db().time().userUnits(), unitDegrees, thetaOuter_());
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
