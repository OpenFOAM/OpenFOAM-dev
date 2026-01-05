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

#include "coneDiskDirectionLagrangianVectorFieldSource.H"
#include "diskInjection.H"
#include "Function1LagrangianFieldSource.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::coneDiskDirectionLagrangianVectorFieldSource::
coneDiskDirectionLagrangianVectorFieldSource
(
    const LagrangianFieldSourceBase& field,
    const dictionary& dict
)
:
    field_(field),
    thetaInner_
    (
        Function1<scalar>::New
        (
            "thetaInner",
            field.db().time().userUnits(),
            unitDegrees,
            dict
        )
    ),
    thetaOuter_
    (
        Function1<scalar>::New
        (
            "thetaOuter",
            field.db().time().userUnits(),
            unitDegrees,
            dict
        )
    )
{}


Foam::coneDiskDirectionLagrangianVectorFieldSource::
coneDiskDirectionLagrangianVectorFieldSource
(
    const coneDiskDirectionLagrangianVectorFieldSource& cddlvfs,
    const LagrangianFieldSourceBase& field
)
:
    field_(field),
    thetaInner_(cddlvfs.thetaInner_, false),
    thetaOuter_(cddlvfs.thetaOuter_, false)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::coneDiskDirectionLagrangianVectorFieldSource::
~coneDiskDirectionLagrangianVectorFieldSource()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::LagrangianSubVectorField>
Foam::coneDiskDirectionLagrangianVectorFieldSource::direction
(
    const LagrangianInjection& injection,
    const LagrangianSubMesh& subMesh
) const
{
    const Lagrangian::diskInjection& diskInjection =
        field_.modelCast<Lagrangian::diskInjection>(injection);

    // Evaluate the cone angle
    const tmp<LagrangianSubScalarField> tthetaInner =
        Function1LagrangianFieldSource::value(subMesh, dimless, thetaInner_());
    const tmp<LagrangianSubScalarField> tthetaOuter =
        Function1LagrangianFieldSource::value(subMesh, dimless, thetaOuter_());
    const LagrangianSubScalarField theta
    (
        (1 - diskInjection.rFrac())*tthetaInner
      + diskInjection.rFrac()*tthetaOuter
    );

    // Return the direction
    return cos(theta)*diskInjection.axis() + sin(theta)*diskInjection.radial();
}


void Foam::coneDiskDirectionLagrangianVectorFieldSource::write
(
    Ostream& os
) const
{
    writeEntry(os, field_.db().time().userUnits(), unitDegrees, thetaInner_());
    writeEntry(os, field_.db().time().userUnits(), unitDegrees, thetaOuter_());
}


// ************************************************************************* //
