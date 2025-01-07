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

#include "coneVelocityLagrangianVectorFieldSource.H"
#include "mathematicalConstants.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::coneVelocityLagrangianVectorFieldSource::
coneVelocityLagrangianVectorFieldSource
(
    const regIOobject& iIo,
    const dictionary& dict
)
:
    LagrangianVectorFieldSource(iIo, dict),
    Function1LagrangianFieldSource<vector>(*this),
    Umean_
    (
        Function1<vector>::New
        (
            "Umean",
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
    ),
    rndGen_
    (
        "rndGen",
        dict,
        randomGenerator::seed(iIo.name() + ':' + dict.dictName()),
        false
    ),
    timeIndex_(-1)
{}


Foam::coneVelocityLagrangianVectorFieldSource::
coneVelocityLagrangianVectorFieldSource
(
    const coneVelocityLagrangianVectorFieldSource& field,
    const regIOobject& iIo
)
:
    LagrangianVectorFieldSource(field, iIo),
    Function1LagrangianFieldSource<vector>(*this),
    Umean_(field.Umean_, false),
    thetaInner_(field.thetaInner_, false),
    thetaOuter_(field.thetaOuter_, false),
    rndGen_(field.rndGen_),
    timeIndex_(-1)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::coneVelocityLagrangianVectorFieldSource::
~coneVelocityLagrangianVectorFieldSource()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::LagrangianSubVectorField>
Foam::coneVelocityLagrangianVectorFieldSource::value
(
    const LagrangianInjection& injection,
    const LagrangianSubMesh& subMesh
) const
{
    // Restart the generator if necessary and set the time index up to date
    rndGen_.start(timeIndex_ == db().time().timeIndex());
    timeIndex_ = db().time().timeIndex();

    // Evaluate mean velocity
    const LagrangianSubVectorField Umean(value(injection, subMesh, Umean_()));

    // Normalise to get the cone axis
    const LagrangianSubVectorField nDir(normalised(Umean));

    // Construct a random direction perpendicular to the cone axis
    const tmp<LagrangianSubVectorField> tt1Dir(normalised(perpendicular(nDir)));
    const tmp<LagrangianSubVectorField> tt2Dir(nDir ^ tt1Dir());
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
        value(injection, subMesh, dimless, thetaInner_());
    const tmp<LagrangianSubScalarField> tthetaOuter =
        value(injection, subMesh, dimless, thetaOuter_());
    const tmp<LagrangianSubScalarField> tthetaFrac =
        LagrangianSubScalarField::New
        (
            "thetaFrac",
            subMesh,
            dimless,
            rndGen_.scalar01(subMesh.size())
        );
    const LagrangianSubScalarField theta
    (
        sqrt
        (
            (1 - tthetaFrac())*sqr(tthetaInner)
          + tthetaFrac()*sqr(tthetaOuter)
        )
    );
    tthetaFrac.clear();

    // Return the velocity in the calculated direction
    return mag(Umean)*(cos(theta)*nDir + sin(theta)*tDir);
}


void Foam::coneVelocityLagrangianVectorFieldSource::write(Ostream& os) const
{
    LagrangianVectorFieldSource::write(os);

    writeEntry(os, db().time().userUnits(), dimVelocity, Umean_());
    writeEntry(os, db().time().userUnits(), unitDegrees, thetaInner_());
    writeEntry(os, db().time().userUnits(), unitDegrees, thetaOuter_());

    writeEntry(os, "rndGen", rndGen_);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makeLagrangianTypeFieldSource
    (
        LagrangianVectorFieldSource,
        coneVelocityLagrangianVectorFieldSource
    );
}

// ************************************************************************* //
