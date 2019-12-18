/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
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

#include "waveDisplacementPointPatchVectorField.H"
#include "pointPatchFields.H"
#include "addToRunTimeSelectionTable.H"
#include "Time.H"
#include "polyMesh.H"
#include "OneConstant.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::waveDisplacementPointPatchVectorField::
waveDisplacementPointPatchVectorField
(
    const pointPatch& p,
    const DimensionedField<vector, pointMesh>& iF
)
:
    fixedValuePointPatchField<vector>(p, iF),
    amplitude_(Zero),
    omega_(0.0),
    waveNumber_(Zero)
{}


Foam::waveDisplacementPointPatchVectorField::
waveDisplacementPointPatchVectorField
(
    const pointPatch& p,
    const DimensionedField<vector, pointMesh>& iF,
    const dictionary& dict
)
:
    fixedValuePointPatchField<vector>(p, iF, dict),
    amplitude_(dict.lookup("amplitude")),
    omega_(dict.lookup<scalar>("omega")),
    waveNumber_(dict.lookupOrDefault<vector>("waveNumber", Zero)),
    startRamp_
    (
        dict.found("startRamp")
      ? Function1<scalar>::New("startRamp", dict)
      : autoPtr<Function1<scalar>>
        (
            new Function1s::OneConstant<scalar>("startRamp")
        )
    ),
    endRamp_
    (
        dict.found("endRamp")
      ? Function1<scalar>::New("endRamp", dict)
      : autoPtr<Function1<scalar>>
        (
            new Function1s::OneConstant<scalar>("endRamp")
        )
    ),
    timeRamp_
    (
        dict.found("timeRamp")
      ? Function1<scalar>::New("timeRamp", dict)
      : autoPtr<Function1<scalar>>
        (
            new Function1s::OneConstant<scalar>("timeRamp")
        )
    )
{
    if (!dict.found("value"))
    {
        updateCoeffs();
    }
}


Foam::waveDisplacementPointPatchVectorField::
waveDisplacementPointPatchVectorField
(
    const waveDisplacementPointPatchVectorField& ptf,
    const pointPatch& p,
    const DimensionedField<vector, pointMesh>& iF,
    const pointPatchFieldMapper& mapper
)
:
    fixedValuePointPatchField<vector>(ptf, p, iF, mapper),
    amplitude_(ptf.amplitude_),
    omega_(ptf.omega_),
    waveNumber_(ptf.waveNumber_)
{}


Foam::waveDisplacementPointPatchVectorField::
waveDisplacementPointPatchVectorField
(
    const waveDisplacementPointPatchVectorField& ptf,
    const DimensionedField<vector, pointMesh>& iF
)
:
    fixedValuePointPatchField<vector>(ptf, iF),
    amplitude_(ptf.amplitude_),
    omega_(ptf.omega_),
    waveNumber_(ptf.waveNumber_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::waveDisplacementPointPatchVectorField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    const polyMesh& mesh = this->internalField().mesh()();
    const Time& t = mesh.time();

    const scalarField points(waveNumber_ & patch().localPoints());

    const scalar timeRamp = timeRamp_->value(t.value());

    const scalarField startRamp(startRamp_->value(points));

    const scalarField endRamp
    (
        endRamp_->value(points[points.size() - 1] - points)
    );

    Field<vector>::operator=
    (
        timeRamp*startRamp*endRamp*amplitude_*cos(omega_*t.value() - points)
    );

    fixedValuePointPatchField<vector>::updateCoeffs();
}


void Foam::waveDisplacementPointPatchVectorField::write(Ostream& os) const
{
    pointPatchField<vector>::write(os);
    writeEntry(os, "amplitude", amplitude_);
    writeEntry(os, "omega", omega_);
    writeEntry(os, "waveNumber", waveNumber_);

    if (!isType<Function1s::OneConstant<scalar>>(startRamp_()))
    {
        writeEntry(os, startRamp_());
    }

    if (!isType<Function1s::OneConstant<scalar>>(endRamp_()))
    {
        writeEntry(os, endRamp_());
    }

    if (!isType<Function1s::OneConstant<scalar>>(timeRamp_()))
    {
        writeEntry(os, timeRamp_());
    }

    writeEntry(os, "value", *this);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePointPatchTypeField
    (
        pointPatchVectorField,
        waveDisplacementPointPatchVectorField
    );
}

// ************************************************************************* //
