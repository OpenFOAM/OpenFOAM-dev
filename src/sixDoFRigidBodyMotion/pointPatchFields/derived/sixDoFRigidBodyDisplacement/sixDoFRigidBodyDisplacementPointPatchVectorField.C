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

#include "sixDoFRigidBodyDisplacementPointPatchVectorField.H"
#include "pointPatchFields.H"
#include "addToRunTimeSelectionTable.H"
#include "Time.H"
#include "fvMesh.H"
#include "volFields.H"
#include "uniformDimensionedFields.H"
#include "forces.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

sixDoFRigidBodyDisplacementPointPatchVectorField::
sixDoFRigidBodyDisplacementPointPatchVectorField
(
    const pointPatch& p,
    const DimensionedField<vector, pointMesh>& iF
)
:
    fixedValuePointPatchField<vector>(p, iF),
    motion_(),
    initialPoints_(p.localPoints()),
    rhoInf_(1.0),
    rhoName_("rho"),
    lookupGravity_(-1),
    g_(vector::zero),
    curTimeIndex_(-1)
{}


sixDoFRigidBodyDisplacementPointPatchVectorField::
sixDoFRigidBodyDisplacementPointPatchVectorField
(
    const pointPatch& p,
    const DimensionedField<vector, pointMesh>& iF,
    const dictionary& dict
)
:
    fixedValuePointPatchField<vector>(p, iF, dict),
    motion_(dict, dict),
    rhoInf_(1.0),
    rhoName_(dict.lookupOrDefault<word>("rhoName", "rho")),
    lookupGravity_(-1),
    g_(vector::zero),
    curTimeIndex_(-1)
{
    if (rhoName_ == "rhoInf")
    {
        rhoInf_ = readScalar(dict.lookup("rhoInf"));
    }

    if (dict.readIfPresent("g", g_))
    {
        lookupGravity_ = -2;
    }

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


sixDoFRigidBodyDisplacementPointPatchVectorField::
sixDoFRigidBodyDisplacementPointPatchVectorField
(
    const sixDoFRigidBodyDisplacementPointPatchVectorField& ptf,
    const pointPatch& p,
    const DimensionedField<vector, pointMesh>& iF,
    const pointPatchFieldMapper& mapper
)
:
    fixedValuePointPatchField<vector>(ptf, p, iF, mapper),
    motion_(ptf.motion_),
    initialPoints_(ptf.initialPoints_, mapper),
    rhoInf_(ptf.rhoInf_),
    rhoName_(ptf.rhoName_),
    lookupGravity_(ptf.lookupGravity_),
    g_(ptf.g_),
    curTimeIndex_(-1)
{}


sixDoFRigidBodyDisplacementPointPatchVectorField::
sixDoFRigidBodyDisplacementPointPatchVectorField
(
    const sixDoFRigidBodyDisplacementPointPatchVectorField& ptf,
    const DimensionedField<vector, pointMesh>& iF
)
:
    fixedValuePointPatchField<vector>(ptf, iF),
    motion_(ptf.motion_),
    initialPoints_(ptf.initialPoints_),
    rhoInf_(ptf.rhoInf_),
    rhoName_(ptf.rhoName_),
    lookupGravity_(ptf.lookupGravity_),
    g_(ptf.g_),
    curTimeIndex_(-1)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void sixDoFRigidBodyDisplacementPointPatchVectorField::autoMap
(
    const pointPatchFieldMapper& m
)
{
    fixedValuePointPatchField<vector>::autoMap(m);

    initialPoints_.autoMap(m);
}


void sixDoFRigidBodyDisplacementPointPatchVectorField::rmap
(
    const pointPatchField<vector>& ptf,
    const labelList& addr
)
{
    const sixDoFRigidBodyDisplacementPointPatchVectorField& sDoFptf =
        refCast<const sixDoFRigidBodyDisplacementPointPatchVectorField>(ptf);

    fixedValuePointPatchField<vector>::rmap(sDoFptf, addr);

    initialPoints_.rmap(sDoFptf.initialPoints_, addr);
}


void sixDoFRigidBodyDisplacementPointPatchVectorField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    if (lookupGravity_ < 0)
    {
        if (db().foundObject<uniformDimensionedVectorField>("g"))
        {
            if (lookupGravity_ == -2)
            {
                FatalErrorIn
                (
                    "void sixDoFRigidBodyDisplacementPointPatchVectorField"
                    "::updateCoeffs()"
                )
                    << "Specifying the value of g in this boundary condition "
                    << "when g is available from the database is considered "
                    << "a fatal error to avoid the possibility of inconsistency"
                    << exit(FatalError);
            }
            else
            {
                lookupGravity_ = 1;
            }
        }
        else
        {
            lookupGravity_ = 0;
        }
    }

    const polyMesh& mesh = this->dimensionedInternalField().mesh()();
    const Time& t = mesh.time();
    const pointPatch& ptPatch = this->patch();

    // Store the motion state at the beginning of the time-step
    if (curTimeIndex_ != t.timeIndex())
    {
        motion_.newTime();
        curTimeIndex_ = t.timeIndex();
    }

    // Patch force data is valid for the current positions, so
    // calculate the forces on the motion object from this data, then
    // update the positions

    motion_.updatePosition(t.deltaTValue(), t.deltaT0Value());

    dictionary forcesDict;

    forcesDict.add("type", forces::typeName);
    forcesDict.add("patches", wordList(1, ptPatch.name()));
    forcesDict.add("rhoInf", rhoInf_);
    forcesDict.add("rhoName", rhoName_);
    forcesDict.add("CofR", motion_.centreOfRotation());

    forces f("forces", db(), forcesDict);

    f.calcForcesMoment();

    // Get the forces on the patch faces at the current positions

    if (lookupGravity_ == 1)
    {
        uniformDimensionedVectorField g =
            db().lookupObject<uniformDimensionedVectorField>("g");

        g_ = g.value();
    }

    // scalar ramp = min(max((t.value() - 5)/10, 0), 1);
    scalar ramp = 1.0;

    motion_.updateAcceleration
    (
        ramp*(f.forceEff() + motion_.mass()*g_),
        ramp*(f.momentEff() + motion_.mass()*(motion_.momentArm() ^ g_)),
        t.deltaTValue()
    );

    Field<vector>::operator=
    (
        motion_.transform(initialPoints_) - initialPoints_
    );

    fixedValuePointPatchField<vector>::updateCoeffs();
}


void sixDoFRigidBodyDisplacementPointPatchVectorField::write(Ostream& os) const
{
    pointPatchField<vector>::write(os);

    os.writeKeyword("rhoName") << rhoName_ << token::END_STATEMENT << nl;

    if (rhoName_ == "rhoInf")
    {
        os.writeKeyword("rhoInf") << rhoInf_ << token::END_STATEMENT << nl;
    }

    if (lookupGravity_ == 0 || lookupGravity_ == -2)
    {
        os.writeKeyword("g") << g_ << token::END_STATEMENT << nl;
    }

    motion_.write(os);

    initialPoints_.writeEntry("initialPoints", os);

    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePointPatchTypeField
(
    pointPatchVectorField,
    sixDoFRigidBodyDisplacementPointPatchVectorField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
