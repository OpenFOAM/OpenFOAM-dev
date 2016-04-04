/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016 OpenFOAM Foundation
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

#include "rigidBodyModel.H"
#include "masslessBody.H"
#include "nullJoint.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace RBD
{
    defineTypeNameAndDebug(rigidBodyModel, 0);
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::RBD::rigidBodyModel::initializeRootBody()
{
    bodies_.append(new masslessBody);
    lambda_.append(0);
    bodyIDs_.insert("root", 0);
    joints_.append(new joints::null(*this));
    XT_.append(spatialTransform());

    nDoF_ = 0;
    nw_ = 0;

    resizeState();
}


void Foam::RBD::rigidBodyModel::resizeState()
{
    Xlambda_.append(spatialTransform());
    X0_.append(spatialTransform());

    v_.append(Zero);
    a_.append(Zero);
    c_.append(Zero);

    IA_.append(spatialTensor::I);
    pA_.append(Zero);

    S_.append(Zero);
    S1_.append(Zero);
    U_.append(Zero);
    U1_.append(Zero);
    Dinv_.append(Zero);
    u_.append(Zero);
}


// * * * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * * //

Foam::RBD::rigidBodyModel::rigidBodyModel()
:
    g_(Zero)
{
    initializeRootBody();
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::label Foam::RBD::rigidBodyModel::join
(
    const label parentID,
    const spatialTransform& XT,
    const autoPtr<joint>& jointPtr,
    const autoPtr<rigidBody>& bodyPtr
)
{
    // Append the body
    const rigidBody& body = bodyPtr();
    bodies_.append(bodyPtr);
    const label bodyID = nBodies()-1;
    bodyIDs_.insert(body.name(), bodyID);

    // If the parentID refers to a merged body find the parent into which it has
    // been merged and set lambda and XT accordingly
    if (merged(parentID))
    {
        const subBody& sBody = mergedBody(parentID);
        lambda_.append(sBody.parentID());
        XT_.append(XT & sBody.parentXT());
    }
    else
    {
        lambda_.append(parentID);
        XT_.append(XT);
    }

    // Append the joint
    const joint& prevJoint = joints_[joints_.size() - 1];
    joints_.append(jointPtr);
    joint& curJoint = joints_[joints_.size() - 1];
    curJoint.index() = joints_.size() - 1;
    curJoint.qIndex() = prevJoint.qIndex() + prevJoint.nDoF();
    curJoint.wIndex() = prevJoint.wIndex() + prevJoint.nw();

    // Increment the degrees of freedom
    nDoF_ += curJoint.nDoF();
    nw_ += curJoint.nw();

    resizeState();

    return bodyID;
}


Foam::label Foam::RBD::rigidBodyModel::join
(
    const label parentID,
    const spatialTransform& XT,
    const PtrList<joint>& compositeJoint,
    const autoPtr<rigidBody>& bodyPtr
)
{
    label parent = parentID;

    // For all but the final joint in the set add a masslessBody with the
    // joint and transform
    for (label j=0; j<compositeJoint.size()-1; j++)
    {
        parent = join
        (
            parent,
            j == 0 ? XT : spatialTransform(),
            compositeJoint[j].clone(),
            autoPtr<rigidBody>(new masslessBody)
        );
    }

    // For the final joint in the set add the read body
    return join
    (
        parent,
        compositeJoint.size() == 1 ? XT : spatialTransform(),
        compositeJoint[compositeJoint.size()-1].clone(),
        bodyPtr
    );
}


Foam::label Foam::RBD::rigidBodyModel::merge
(
    const label parentID,
    const spatialTransform& XT,
    const autoPtr<rigidBody>& bodyPtr
)
{
    autoPtr<subBody> sBodyPtr;

    // If the parentID refers to a merged body find the parent into which it has
    // been merged and merge this on into the same parent with the appropriate
    // transform
    if (merged(parentID))
    {
        const subBody& sBody = mergedBody(parentID);
        sBodyPtr.set
        (
            new subBody
            (
                bodyPtr,
                bodies_[sBody.parentID()].name(),
                sBody.parentID(),
                XT & sBody.parentXT()
            )
        );
    }
    else
    {
        sBodyPtr.set
        (
            new subBody
            (
                bodyPtr,
                bodies_[parentID].name(),
                parentID,
                XT
            )
        );
    }

    const subBody& sBody = sBodyPtr();
    mergedBodies_.append(sBodyPtr);

    // Merge the sub-body with the parent
    bodies_[sBody.parentID()].merge(sBody);

    const label sBodyID = mergedBodyID(mergedBodies_.size() - 1);
    bodyIDs_.insert(sBody.name(), sBodyID);

    return sBodyID;
}


Foam::spatialTransform Foam::RBD::rigidBodyModel::X0
(
    const label bodyId
) const
{
    if (merged(bodyId))
    {
        const subBody& mBody = mergedBody(bodyId);
        return mBody.parentXT() & X0_[mBody.parentID()];
    }
    else
    {
        return X0_[bodyId];
    }
}


// ************************************************************************* //
