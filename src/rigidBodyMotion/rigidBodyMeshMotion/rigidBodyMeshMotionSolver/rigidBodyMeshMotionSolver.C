/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2025 OpenFOAM Foundation
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

#include "rigidBodyMeshMotionSolver.H"
#include "polyMesh.H"
#include "polyTopoChangeMap.H"
#include "pointConstraints.H"
#include "timeIOdictionary.H"
#include "uniformDimensionedFields.H"
#include "forces.H"
#include "OneConstant.H"
#include "mathematicalConstants.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(rigidBodyMeshMotionSolver, 0);

    addToRunTimeSelectionTable
    (
        motionSolver,
        rigidBodyMeshMotionSolver,
        dictionary
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::rigidBodyMeshMotionSolver::bodyMesh::bodyMesh
(
    const polyMesh& mesh,
    const word& name,
    const label bodyID,
    const dictionary& dict
)
:
    name_(name),
    bodyIndex_(bodyID),
    patches_(wordReList(dict.lookup("patches"))),
    patchSet_(mesh.boundaryMesh().patchSet(patches_))
{}


Foam::rigidBodyMeshMotionSolver::rigidBodyMeshMotionSolver
(
    const word& name,
    const polyMesh& mesh,
    const dictionary& dict
)
:
    motionSolver(name, mesh, typeName),
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
    curTimeIndex_(-1),
    meshSolverPtr_
    (
        motionSolver::New
        (
            name,
            mesh,
            IOdictionary
            (
                IOobject
                (
                    typedName("meshSolver"),
                    mesh.time().constant(),
                    mesh
                ),
                dict.subDict("meshSolver")
            )
        )
    ),
    meshSolver_(refCast<displacementMotionSolver>(meshSolverPtr_()))
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

    const dictionary& bodiesDict = dict.subDict("bodies");

    forAllConstIter(IDLList<entry>, bodiesDict, iter)
    {
        const dictionary& bodyDict = iter().dict();

        if (bodyDict.found("patches"))
        {
            const label bodyID = this->bodyIndex(iter().keyword());

            if (bodyID == -1)
            {
                FatalErrorInFunction
                    << "Body " << iter().keyword()
                    << " has been merged with another body"
                       " and cannot be assigned a set of patches"
                    << exit(FatalError);
            }

            bodyMeshes_.append
            (
                new bodyMesh
                (
                    mesh,
                    iter().keyword(),
                    bodyID,
                    bodyDict
                )
            );
        }
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::rigidBodyMeshMotionSolver::~rigidBodyMeshMotionSolver()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::pointField>
Foam::rigidBodyMeshMotionSolver::curPoints() const
{
    return meshSolverPtr_->curPoints();
}


void Foam::rigidBodyMeshMotionSolver::solve()
{
    const Time& t = mesh().time();

    if (mesh().nPoints() != meshSolver_.points0().size())
    {
        FatalErrorInFunction
            << "The number of points in the mesh seems to have changed." << endl
            << "In constant/polyMesh there are " << meshSolver_.points0().size()
            << " points; in the current mesh there are " << mesh().nPoints()
            << " points." << exit(FatalError);
    }

    // Store the motion state at the beginning of the time-step
    if (curTimeIndex_ != t.timeIndex())
    {
        newTime();
        curTimeIndex_ = t.timeIndex();
    }

    const scalar ramp = ramp_->value(t.value());

    if (mesh().foundObject<uniformDimensionedVectorField>("g"))
    {
        g() =
            ramp
           *mesh().lookupObject<uniformDimensionedVectorField>("g").value();
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
            const label bodyID = bodyMeshes_[bi].bodyIndex_;

            functionObjects::forces f
            (
                functionObjects::forces::typeName,
                t,
                dictionary::entries
                (
                    "type", functionObjects::forces::typeName,
                    "patches", bodyMeshes_[bi].patches_,
                    "rhoInf", rhoInf_,
                    "rho", rhoName_,
                    "CofR", vector::zero
                )
            );

            f.calcForcesMoments();

            fx[bodyID] = ramp*spatialVector(f.momentEff(), f.forceEff());
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
            status(bodyMeshes_[bi].bodyIndex_);
        }
    }

    // Update the displacements
    forAll(bodyMeshes_, bi)
    {
        forAllConstIter(labelHashSet, bodyMeshes_[bi].patchSet_, iter)
        {
            const label patchi = iter.key();

            const pointField patchPoints0
            (
                meshSolver_.pointDisplacement().boundaryField()[patchi]
               .patchInternalField(meshSolver_.points0())
            );

            meshSolver_.pointDisplacement().boundaryFieldRef()[patchi] ==
            (
                Foam::transformPoints
                (
                    transform0(bodyMeshes_[bi].bodyIndex_),
                    patchPoints0
                ) - patchPoints0
            )();
        }
    }

    meshSolverPtr_->solve();
}


void Foam::rigidBodyMeshMotionSolver::movePoints(const pointField& points)
{
    meshSolverPtr_->movePoints(points);
}


void Foam::rigidBodyMeshMotionSolver::topoChange(const polyTopoChangeMap& map)
{
    meshSolverPtr_->topoChange(map);
}


void Foam::rigidBodyMeshMotionSolver::mapMesh(const polyMeshMap& map)
{
    meshSolverPtr_->mapMesh(map);
}


void Foam::rigidBodyMeshMotionSolver::distribute
(
    const polyDistributionMap& map
)
{
    meshSolverPtr_->distribute(map);
}


bool Foam::rigidBodyMeshMotionSolver::write() const
{
    timeIOdictionary dict
    (
        IOobject
        (
            "rigidBodyMotionState",
            mesh().time().name(),
            "uniform",
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        )
    );

    state().write(dict);

    return
        motionSolver::write()
     && dict.regIOobject::writeObject
        (
            IOstream::ASCII,
            IOstream::currentVersion,
            mesh().time().writeCompression(),
            true
        );
}


// ************************************************************************* //
