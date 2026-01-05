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

#include "coneDirectionLagrangianVectorFieldSource.H"
#include "Function1LagrangianFieldSource.H"
#include "mathematicalConstants.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::coneDirectionLagrangianVectorFieldSource::
coneDirectionLagrangianVectorFieldSource
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
    ),
    rndGen_
    (
        "coneDirectionRndGen",
        dict,
        randomGenerator::seed(field.db().name() + ':' + dict.dictName()),
        false
    ),
    timeIndex_(-1)
{}


Foam::coneDirectionLagrangianVectorFieldSource::
coneDirectionLagrangianVectorFieldSource
(
    const coneDirectionLagrangianVectorFieldSource& cdlvfs,
    const LagrangianFieldSourceBase& field
)
:
    field_(field),
    thetaInner_(cdlvfs.thetaInner_, false),
    thetaOuter_(cdlvfs.thetaOuter_, false),
    rndGen_(cdlvfs.rndGen_),
    timeIndex_(-1)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::coneDirectionLagrangianVectorFieldSource::
~coneDirectionLagrangianVectorFieldSource()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::LagrangianSubVectorField>
Foam::coneDirectionLagrangianVectorFieldSource::direction
(
    const LagrangianSubVectorField& axis
) const
{
    const LagrangianSubMesh& subMesh = axis.mesh();

    // Restart the generator if necessary and set the time index up to date
    rndGen_.start(timeIndex_ == field_.db().time().timeIndex());
    timeIndex_ = field_.db().time().timeIndex();

    // Construct a random direction perpendicular to the cone axis
    const tmp<LagrangianSubVectorField> tt1Dir(normalised(perpendicular(axis)));
    const tmp<LagrangianSubVectorField> tt2Dir(axis ^ tt1Dir());
    const tmp<LagrangianSubScalarField> tphi =
        LagrangianSubScalarField::New
        (
            "phi",
            subMesh,
            dimless,
            constant::mathematical::twoPi*rndGen_.scalar01(subMesh.size())
        );
    const LagrangianSubVectorField tDir
    (
        tt1Dir*cos(tphi()) + tt2Dir*sin(tphi())
    );
    tphi.clear();

    // Pick a random angle within the cone angles
    const tmp<LagrangianSubScalarField> tthetaInner =
        Function1LagrangianFieldSource::value(subMesh, dimless, thetaInner_());
    const tmp<LagrangianSubScalarField> tthetaOuter =
        Function1LagrangianFieldSource::value(subMesh, dimless, thetaOuter_());
    const tmp<LagrangianSubScalarField> tfrac =
        LagrangianSubScalarField::New
        (
            "frac",
            subMesh,
            dimless,
            rndGen_.scalar01(subMesh.size())
        );
    const LagrangianSubScalarField theta
    (
        sqrt
        (
            (1 - tfrac())*sqr(tthetaInner)
          + tfrac()*sqr(tthetaOuter)
        )
    );
    tfrac.clear();

    // Return the direction
    return cos(theta)*axis + sin(theta)*tDir;
}


void Foam::coneDirectionLagrangianVectorFieldSource::write(Ostream& os) const
{
    writeEntry(os, field_.db().time().userUnits(), unitDegrees, thetaInner_());
    writeEntry(os, field_.db().time().userUnits(), unitDegrees, thetaOuter_());

    writeEntry(os, "coneDirectionRndGen", rndGen_);
}


// ************************************************************************* //
