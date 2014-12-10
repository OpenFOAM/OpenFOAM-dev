/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2014 OpenFOAM Foundation
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

#include "uncoupledSixDoFRigidBodyDisplacementPointPatchVectorField.H"
#include "pointPatchFields.H"
#include "addToRunTimeSelectionTable.H"
#include "Time.H"
#include "fvMesh.H"
#include "volFields.H"
#include "uniformDimensionedFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

uncoupledSixDoFRigidBodyDisplacementPointPatchVectorField::
uncoupledSixDoFRigidBodyDisplacementPointPatchVectorField
(
    const pointPatch& p,
    const DimensionedField<vector, pointMesh>& iF
)
:
    fixedValuePointPatchField<vector>(p, iF),
    motion_(),
    initialPoints_(p.localPoints()),
    curTimeIndex_(-1)
{}


uncoupledSixDoFRigidBodyDisplacementPointPatchVectorField::
uncoupledSixDoFRigidBodyDisplacementPointPatchVectorField
(
    const pointPatch& p,
    const DimensionedField<vector, pointMesh>& iF,
    const dictionary& dict
)
:
    fixedValuePointPatchField<vector>(p, iF, dict),
    motion_(dict, dict),
    curTimeIndex_(-1)
{
    if (!dict.found("value"))
    {
        updateCoeffs();
    }

    if (dict.found("initialPoints"))
    {
        initialPoints_ = vectorField("initialPoints", dict , p.size());
    }
    else
    {
        initialPoints_ = p.localPoints();
    }
}


uncoupledSixDoFRigidBodyDisplacementPointPatchVectorField::
uncoupledSixDoFRigidBodyDisplacementPointPatchVectorField
(
    const uncoupledSixDoFRigidBodyDisplacementPointPatchVectorField& ptf,
    const pointPatch& p,
    const DimensionedField<vector, pointMesh>& iF,
    const pointPatchFieldMapper& mapper
)
:
    fixedValuePointPatchField<vector>(ptf, p, iF, mapper),
    motion_(ptf.motion_),
    initialPoints_(ptf.initialPoints_, mapper),
    curTimeIndex_(-1)
{}


uncoupledSixDoFRigidBodyDisplacementPointPatchVectorField::
uncoupledSixDoFRigidBodyDisplacementPointPatchVectorField
(
    const uncoupledSixDoFRigidBodyDisplacementPointPatchVectorField& ptf,
    const DimensionedField<vector, pointMesh>& iF
)
:
    fixedValuePointPatchField<vector>(ptf, iF),
    motion_(ptf.motion_),
    initialPoints_(ptf.initialPoints_),
    curTimeIndex_(-1)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void uncoupledSixDoFRigidBodyDisplacementPointPatchVectorField::autoMap
(
    const pointPatchFieldMapper& m
)
{
    fixedValuePointPatchField<vector>::autoMap(m);

    initialPoints_.autoMap(m);
}


void uncoupledSixDoFRigidBodyDisplacementPointPatchVectorField::rmap
(
    const pointPatchField<vector>& ptf,
    const labelList& addr
)
{
    const uncoupledSixDoFRigidBodyDisplacementPointPatchVectorField& uSDoFptf =
    refCast
    <
        const uncoupledSixDoFRigidBodyDisplacementPointPatchVectorField
    >(ptf);

    fixedValuePointPatchField<vector>::rmap(uSDoFptf, addr);

    initialPoints_.rmap(uSDoFptf.initialPoints_, addr);
}


void uncoupledSixDoFRigidBodyDisplacementPointPatchVectorField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    const polyMesh& mesh = this->dimensionedInternalField().mesh()();
    const Time& t = mesh.time();

    // Store the motion state at the beginning of the time-step
    if (curTimeIndex_ != t.timeIndex())
    {
        motion_.newTime();
        curTimeIndex_ = t.timeIndex();
    }

    motion_.updatePosition(t.deltaTValue(), t.deltaT0Value());

    vector gravity = vector::zero;

    if (db().foundObject<uniformDimensionedVectorField>("g"))
    {
        uniformDimensionedVectorField g =
        db().lookupObject<uniformDimensionedVectorField>("g");

        gravity = g.value();
    }

    // Do not modify the accelerations, except with gravity, where available
    motion_.updateAcceleration
    (
        gravity*motion_.mass(),
        vector::zero,
        t.deltaTValue()
    );

    Field<vector>::operator=
    (
        motion_.transform(initialPoints_) - initialPoints_
    );

    fixedValuePointPatchField<vector>::updateCoeffs();
}


void uncoupledSixDoFRigidBodyDisplacementPointPatchVectorField::write
(
    Ostream& os
) const
{
    pointPatchField<vector>::write(os);
    motion_.write(os);
    initialPoints_.writeEntry("initialPoints", os);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePointPatchTypeField
(
    pointPatchVectorField,
    uncoupledSixDoFRigidBodyDisplacementPointPatchVectorField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
