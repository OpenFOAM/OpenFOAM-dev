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

#include "fanDirectionLagrangianVectorFieldSource.H"
#include "Function1LagrangianFieldSource.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fanDirectionLagrangianVectorFieldSource::
fanDirectionLagrangianVectorFieldSource
(
    const LagrangianFieldSourceBase& field,
    const dictionary& dict
)
:
    field_(field),
    normal_
    (
        Function1<vector>::New
        (
            "normal",
            field.db().time().userUnits(),
            dimless,
            dict
        )
    ),
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
    ),
    rndGen_
    (
        "fanDirectionRndGen",
        dict,
        randomGenerator::seed(field.db().name() + ':' + dict.dictName()),
        false
    ),
    timeIndex_(-1)
{}


Foam::fanDirectionLagrangianVectorFieldSource::
fanDirectionLagrangianVectorFieldSource
(
    const fanDirectionLagrangianVectorFieldSource& cdlvfs,
    const LagrangianFieldSourceBase& field
)
:
    field_(field),
    normal_(cdlvfs.normal_, false),
    thetaInner_(cdlvfs.thetaInner_, false),
    thetaOuter_(cdlvfs.thetaOuter_, false),
    rndGen_(cdlvfs.rndGen_),
    timeIndex_(-1)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fanDirectionLagrangianVectorFieldSource::
~fanDirectionLagrangianVectorFieldSource()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::LagrangianSubVectorField>
Foam::fanDirectionLagrangianVectorFieldSource::direction
(
    const LagrangianInjection& injection,
    const LagrangianSubVectorField& axis
) const
{
    const LagrangianSubMesh& subMesh = axis.mesh();

    // Restart the generator if necessary and set the time index up to date
    rndGen_.start(timeIndex_ == field_.db().time().timeIndex());
    timeIndex_ = field_.db().time().timeIndex();

    // Construct a direction in the plane of the fan
    const tmp<LagrangianSubVectorField> normal =
        Function1LagrangianFieldSource::value
        (
            injection,
            subMesh,
            dimless,
            normal_()
        );
    const tmp<LagrangianSubVectorField> tangential(axis ^ normal);

    // Pick a random angle within the fan angles
    const tmp<LagrangianSubScalarField> tthetaInner =
        Function1LagrangianFieldSource::value
        (
            injection,
            subMesh,
            dimless,
            thetaInner_()
        );
    const tmp<LagrangianSubScalarField> tthetaOuter =
        Function1LagrangianFieldSource::value
        (
            injection,
            subMesh,
            dimless,
            thetaOuter_()
        );
    const tmp<LagrangianSubScalarField> tfrac =
        LagrangianSubScalarField::New
        (
            "frac",
            subMesh,
            dimless,
            rndGen_.scalarAB(subMesh.size(), -1, 1)
        );
    const tmp<LagrangianSubScalarField> tmagFrac = mag(tfrac());
    const LagrangianSubScalarField theta
    (
        (1 - tmagFrac())*tthetaInner
      + tmagFrac()*tthetaOuter
    );
    tmagFrac.clear();

    // Return the calculated direction
    return cos(theta)*axis + sign(tfrac)*sin(theta)*tangential;
}


void Foam::fanDirectionLagrangianVectorFieldSource::write(Ostream& os) const
{
    writeEntry(os, field_.db().time().userUnits(), unitNone, normal_());
    writeEntry(os, field_.db().time().userUnits(), unitDegrees, thetaInner_());
    writeEntry(os, field_.db().time().userUnits(), unitDegrees, thetaOuter_());

    writeEntry(os, "fanDirectionRndGen", rndGen_);
}


// ************************************************************************* //
