/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2019 OpenFOAM Foundation
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
#include "addToRunTimeSelectionTable.H"
#include "polyMesh.H"
#include "pointPatchDist.H"
#include "pointConstraints.H"
#include "uniformDimensionedFields.H"
#include "forces.H"
#include "OneConstant.H"
#include "mathematicalConstants.H"

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
    bodyID_(bodyID),
    patches_(wordReList(dict.lookup("patches"))),
    patchSet_(mesh.boundaryMesh().patchSet(patches_))
{}


Foam::rigidBodyMeshMotionSolver::rigidBodyMeshMotionSolver
(
    const polyMesh& mesh,
    const dictionary& dict
)
:
    motionSolver(mesh, dict, typeName),
    RBD::rigidBodyMotion
    (
        coeffDict(),
        IOobject
        (
            "rigidBodyMotionState",
            mesh.time().timeName(),
            "uniform",
            mesh
        ).typeHeaderOk<IOdictionary>(true)
      ? IOdictionary
        (
            IOobject
            (
                "rigidBodyMotionState",
                mesh.time().timeName(),
                "uniform",
                mesh,
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE,
                false
            )
        )
      : coeffDict()
    ),
    test_(coeffDict().lookupOrDefault<Switch>("test", false)),
    rhoInf_(1.0),
    rhoName_(coeffDict().lookupOrDefault<word>("rho", "rho")),
    ramp_(nullptr),
    curTimeIndex_(-1),
    meshSolverPtr_
    (
        motionSolver::New
        (
            mesh,
            IOdictionary
            (
                IOobject
                (
                    "rigidBodyMotionSolver:meshSolver",
                    mesh.time().constant(),
                    mesh
                ),
                coeffDict().subDict("meshSolver")
            )
        )
    ),
    meshSolver_(refCast<displacementMotionSolver>(meshSolverPtr_()))
{
    if (rhoName_ == "rhoInf")
    {
        rhoInf_ = coeffDict().lookup<scalar>("rhoInf");
    }

    if (coeffDict().found("ramp"))
    {
        ramp_ = Function1<scalar>::New("ramp", coeffDict());
    }
    else
    {
        ramp_ = new Function1s::OneConstant<scalar>("ramp");
    }

    const dictionary& bodiesDict = coeffDict().subDict("bodies");

    forAllConstIter(IDLList<entry>, bodiesDict, iter)
    {
        const dictionary& bodyDict = iter().dict();

        if (bodyDict.found("patches"))
        {
            const label bodyID = this->bodyID(iter().keyword());

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
        label nIter(coeffDict().lookup<label>("nIter"));

        for (label i=0; i<nIter; i++)
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
            const label bodyID = bodyMeshes_[bi].bodyID_;

            dictionary forcesDict;
            forcesDict.add("type", functionObjects::forces::typeName);
            forcesDict.add("patches", bodyMeshes_[bi].patches_);
            forcesDict.add("rhoInf", rhoInf_);
            forcesDict.add("rho", rhoName_);
            forcesDict.add("CofR", vector::zero);

            functionObjects::forces f("forces", t, forcesDict);
            f.calcForcesMoment();

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
            status(bodyMeshes_[bi].bodyID_);
        }
    }

    // Update the displacements
    forAll(bodyMeshes_, bi)
    {
        forAllConstIter(labelHashSet, bodyMeshes_[bi].patchSet_, iter)
        {
            label patchi = iter.key();

            pointField patchPoints0
            (
                meshSolver_.pointDisplacement().boundaryField()[patchi]
               .patchInternalField(meshSolver_.points0())
            );

            meshSolver_.pointDisplacement().boundaryFieldRef()[patchi] ==
            (
                transformPoints
                (
                    bodyMeshes_[bi].bodyID_,
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


void Foam::rigidBodyMeshMotionSolver::updateMesh(const mapPolyMesh& mpm)
{
    meshSolverPtr_->updateMesh(mpm);
}


bool Foam::rigidBodyMeshMotionSolver::write() const
{
    IOdictionary dict
    (
        IOobject
        (
            "rigidBodyMotionState",
            mesh().time().timeName(),
            "uniform",
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        )
    );

    state().write(dict);

    return
        dict.regIOobject::writeObject
        (
            IOstream::ASCII,
            IOstream::currentVersion,
            mesh().time().writeCompression(),
            true
        )
     && motionSolver::write();
}


// ************************************************************************* //
