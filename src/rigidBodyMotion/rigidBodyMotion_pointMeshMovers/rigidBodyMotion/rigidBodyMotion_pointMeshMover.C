/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2026 OpenFOAM Foundation
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

#include "rigidBodyMotion_pointMeshMover.H"
#include "timeIOdictionary.H"
#include "uniformDimensionedFields.H"
#include "forces.H"
#include "OneConstant.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace pointMeshMovers
{

    defineTypeNameAndDebug(rigidBodyMotion, 0);

    addToRunTimeSelectionTable
    (
        pointMeshMover,
        rigidBodyMotion,
        dictionary
    );
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::List<Foam::septernion>
Foam::pointMeshMovers::rigidBodyMotion::transforms0() const
{
    List<septernion> transforms0(bodyMeshes_.size());

    forAll(bodyMeshes_, bi)
    {
        if (bodyMeshes_[bi].bodyIndex != -1)
        {
            // Calculate the septernion equivalent of the transformation
            transforms0[bi] =
                septernion(transform0(bodyMeshes_[bi].bodyIndex));
        }
        else
        {
            // Assume the external body is stationary
            transforms0[bi] = septernion::I;
        }
    }

    return transforms0;
}


void Foam::pointMeshMovers::rigidBodyMotion::moveBodies()
{
    const Time& t = poly().time();

    // Store the motion state at the beginning of the time-step
    if (curTimeIndex_ != t.timeIndex())
    {
        newTime();
        curTimeIndex_ = t.timeIndex();
    }

    const scalar ramp = ramp_->value(t.value());

    if (poly().foundObject<uniformDimensionedVectorField>("g"))
    {
        g() =
            ramp
           *poly().lookupObject<uniformDimensionedVectorField>("g").value();
    }

    if (test_)
    {
        for (label i=0; i<nIter_; i++)
        {
            RBD::rigidBodyMotion::solve
            (
                t.value(),
                t.deltaTValue(),
                scalarField(nDoF(), Zero),
                Field<spatialVector>(nBodies(), Zero)
            );
        }
    }
    else
    {
        Field<spatialVector> fx(nBodies(), Zero);

        forAll(bodyMeshes_, bi)
        {
            const label bodyID = bodyMeshes_[bi].bodyIndex;

            if (bodyID != -1)
            {
                functionObjects::forces f
                (
                    functionObjects::forces::typeName,
                    t,
                    dictionary::entries
                    (
                        "type", functionObjects::forces::typeName,
                        "patches", bodyMeshes_[bi].patches(),
                        "rhoInf", rhoInf_,
                        "rho", rhoName_,
                        "CofR", vector::zero
                    )
                );

                f.calcForcesMoments();

                fx[bodyID] = ramp*spatialVector(f.momentEff(), f.forceEff());
            }
        }

        RBD::rigidBodyMotion::solve
        (
            t.value(),
            t.deltaTValue(),
            scalarField(nDoF(), Zero),
            fx
        );
    }

    if (Pstream::master() && report())
    {
        forAll(bodyMeshes_, bi)
        {
            if (bodyMeshes_[bi].bodyIndex != -1)
            {
                status(bodyMeshes_[bi].bodyIndex);
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::pointMeshMovers::rigidBodyMotion::rigidBodyMotion
(
    const polyMesh& mesh,
    const dictionary& dict
)
:
    pointMeshMovers::multiRigidBody(mesh, dict),
    RBD::rigidBodyMotion
    (
        dict,
        typeIOobject<timeIOdictionary>
        (
            "rigidBodyMotionState",
            mesh.time().name(),
            "uniform",
            mesh
        ).headerOk()
      ? timeIOdictionary
        (
            IOobject
            (
                "rigidBodyMotionState",
                mesh.time().name(),
                "uniform",
                mesh,
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE,
                false
            )
        )
      : dict
    ),
    test_(dict.lookupOrDefault<Switch>("test", false)),
    nIter_(test_ ? dict.lookup<label>("nIter") : 0),
    rhoInf_(1.0),
    rhoName_(dict.lookupOrDefault<word>("rho", "rho")),
    ramp_(nullptr),
    curTimeIndex_(-1)
{
    if (rhoName_ == "rhoInf")
    {
        rhoInf_ = dict.lookup<scalar>("rhoInf");
    }

    if (dict.found("ramp"))
    {
        ramp_ = Function1<scalar>::New("ramp", dimTime, dimless, dict);
    }
    else
    {
        ramp_ = new Function1s::OneConstant<scalar>("ramp");
    }

    forAll(bodyMeshes_, bi)
    {
        if (bodyMeshes_[bi].name() != "exterior")
        {
            const label bodyID = this->bodyIndex(bodyMeshes_[bi].name());

            if (bodyID == -1)
            {
                FatalErrorInFunction
                    << "Body " << bodyMeshes_[bi].name()
                    << " has been merged with another body"
                       " and cannot be assigned a set of patches"
                    << exit(FatalError);
            }

            bodyMeshes_[bi].bodyIndex = bodyID;
        }
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::pointMeshMovers::rigidBodyMotion::~rigidBodyMotion()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::pointMeshMovers::rigidBodyMotion::write() const
{
    timeIOdictionary dict
    (
        IOobject
        (
            "rigidBodyMotionState",
            poly().time().name(),
            "uniform",
            poly(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        )
    );

    state().write(dict);

    return
        pointMeshMovers::multiRigidBody::write()
     && dict.regIOobject::writeObject
        (
            IOstream::ASCII,
            IOstream::currentVersion,
            poly().time().writeCompression(),
            true
        );
}


// ************************************************************************* //
