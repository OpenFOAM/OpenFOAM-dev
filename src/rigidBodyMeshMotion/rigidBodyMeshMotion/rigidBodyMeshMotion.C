/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2022 OpenFOAM Foundation
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

#include "rigidBodyMeshMotion.H"
#include "polyMesh.H"
#include "polyTopoChangeMap.H"
#include "pointPatchDist.H"
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
    defineTypeNameAndDebug(rigidBodyMeshMotion, 0);

    addToRunTimeSelectionTable
    (
        motionSolver,
        rigidBodyMeshMotion,
        dictionary
    );
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::List<Foam::septernion> Foam::rigidBodyMeshMotion::transforms0() const
{
    List<septernion> transforms0(bodyMeshes_.size() + 1);

    forAll(bodyMeshes_, bi)
    {
        // Calculate the septernion equivalent of the transformation
        transforms0[bi] = septernion(transform0(bodyMeshes_[bi].bodyID_));
    }

    transforms0[bodyMeshes_.size()] = septernion::I;

    return transforms0;
}


Foam::List<Foam::scalar>& Foam::rigidBodyMeshMotion::weights
(
    const label pointi,
    List<scalar>& w
) const
{
    // Initialise to 1 for the far-field weight
    scalar sum1mw = 1;

    forAll(bodyMeshes_, bi)
    {
        w[bi] = bodyMeshes_[bi].weight_[pointi];
        sum1mw += w[bi]/(1 + small - w[bi]);
    }

    // Calculate the limiter for wi/(1 - wi) to ensure the sum(wi) = 1
    scalar lambda = 1/sum1mw;

    // Limit wi/(1 - wi) and sum the resulting wi
    scalar sumw = 0;
    forAll(bodyMeshes_, bi)
    {
        w[bi] = lambda*w[bi]/(1 + small - w[bi]);
        sumw += w[bi];
    }

    // Calculate the weight for the stationary far-field
    w[bodyMeshes_.size()] = 1 - sumw;

    return w;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::rigidBodyMeshMotion::bodyMesh::bodyMesh
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
    patchSet_(mesh.boundaryMesh().patchSet(patches_)),
    di_(dict.lookup<scalar>("innerDistance")),
    do_(dict.lookup<scalar>("outerDistance")),
    weight_
    (
        IOobject
        (
            name_ + ".motionScale",
            mesh.time().name(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        pointMesh::New(mesh),
        dimensionedScalar(dimless, 0)
    )
{}


Foam::rigidBodyMeshMotion::rigidBodyMeshMotion
(
    const word& name,
    const polyMesh& mesh,
    const dictionary& dict
)
:
    displacementMotionSolver(name, mesh, dict, typeName),
    RBD::rigidBodyMotion
    (
        coeffDict(),
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
      : coeffDict()
    ),
    test_(coeffDict().lookupOrDefault<Switch>("test", false)),
    rhoInf_(1.0),
    rhoName_(coeffDict().lookupOrDefault<word>("rho", "rho")),
    ramp_(nullptr),
    curTimeIndex_(-1)
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

    const pointMesh& pMesh = pointMesh::New(mesh);

    // Calculate scaling factor everywhere for each meshed body
    forAll(bodyMeshes_, bi)
    {
        const pointPatchDist pDist(pMesh, bodyMeshes_[bi].patchSet_, points0());

        bodyMeshes_[bi].weight_.primitiveFieldRef() =
            bodyMeshes_[bi].weight(pDist.primitiveField());

        pointConstraints::New(pMesh).constrain(bodyMeshes_[bi].weight_);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::rigidBodyMeshMotion::~rigidBodyMeshMotion()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Type Foam::rigidBodyMeshMotion::bodyMesh::weight
(
    const Type& pDist
) const
{
    // Scaling: 1 up to di then linear down to 0 at do away from patches
    Type weight
    (
        min(max((do_ - pDist)/(do_ - di_), scalar(0)), scalar(1))
    );

    // Convert the weight function to a cosine
    weight =
        min
        (
            max
            (
                0.5 - 0.5*cos(weight*Foam::constant::mathematical::pi),
                scalar(0)
            ),
            scalar(1)
        );

    return weight;
}


Foam::tmp<Foam::pointField>
Foam::rigidBodyMeshMotion::curPoints() const
{
    return points0() + pointDisplacement_.primitiveField();
}


void Foam::rigidBodyMeshMotion::solve()
{
    const Time& t = mesh().time();

    if (mesh().nPoints() != points0().size())
    {
        FatalErrorInFunction
            << "The number of points in the mesh seems to have changed." << endl
            << "In constant/polyMesh there are " << points0().size()
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

    vectorField& pointDisplacement = pointDisplacement_.primitiveFieldRef();
    const pointField& points0 = this->points0();

    // Update the displacements
    if (bodyMeshes_.size() == 1)
    {
        const septernion transform0(this->transform0(bodyMeshes_[0].bodyID_));
        const scalarField& weight = bodyMeshes_[0].weight_;

        forAll(points0, pointi)
        {
            // Move non-stationary points
            if (weight[pointi] > small)
            {
                // Use solid-body motion where weight = 1
                if (weight[pointi] > 1 - small)
                {
                    pointDisplacement[pointi] =
                        transform0.transformPoint(points0[pointi])
                      - points0[pointi];
                }
                // Slerp septernion interpolation
                else
                {
                    pointDisplacement[pointi] =
                        slerp(septernion::I, transform0, weight[pointi])
                       .transformPoint(points0[pointi])
                      - points0[pointi];
                }
            }
        }
    }
    else
    {
        const List<septernion> transforms0(this->transforms0());
        List<scalar> w(transforms0.size());

        forAll(points0, pointi)
        {
            pointDisplacement[pointi] =
                average(transforms0, weights(pointi, w))
               .transformPoint(points0[pointi])
              - points0[pointi];
        }
    }

    // Displacement has changed. Update boundary conditions
    pointConstraints::New
    (
        pointDisplacement_.mesh()
    ).constrainDisplacement(pointDisplacement_);
}


void Foam::rigidBodyMeshMotion::topoChange(const polyTopoChangeMap& map)
{
    // pointMesh already updates pointFields

    // Get the new points either from the map or the mesh
    const pointField& points =
    (
        map.hasMotionPoints()
      ? map.preMotionPoints()
      : mesh().points()
    );

    const pointMesh& pMesh = pointMesh::New(mesh());
    pointField points0(mesh().points());

    // Iterate to update the transformation of the new points to the
    // corresponding points0, required because the body-point weights are
    // calculated for points0
    for (int iter=0; iter<3; iter++)
    {
        // Calculate scaling factor everywhere for each meshed body
        forAll(bodyMeshes_, bi)
        {
            const pointPatchDist pDist
            (
                pMesh,
                bodyMeshes_[bi].patchSet_,
                points0
            );

            pointScalarField& weight = bodyMeshes_[bi].weight_;

            forAll(points0, pointi)
            {
                const label oldPointi = map.pointMap()[pointi];

                if (oldPointi >= 0)
                {
                    if (map.reversePointMap()[oldPointi] != pointi)
                    {
                        weight[pointi] = bodyMeshes_[bi].weight(pDist[pointi]);
                    }
                }
                else
                {
                    FatalErrorInFunction
                        << "Cannot determine co-ordinates "
                           "of introduced vertices."
                        << " New vertex " << pointi << " at co-ordinate "
                        << points[pointi] << exit(FatalError);
                }
            }

            pointConstraints::New(pMesh).constrain(weight);
        }

        forAll(points0, pointi)
        {
            const label oldPointi = map.pointMap()[pointi];

            if (oldPointi >= 0)
            {
                if (map.reversePointMap()[oldPointi] == pointi)
                {
                    // points0[pointi] = points0_[oldPointi];
                    points0[pointi] = points0_[pointi];
                }
                else
                {
                    if (bodyMeshes_.size() == 1)
                    {
                        // Use solid-body motion where weight = 1
                        if (bodyMeshes_[0].weight_[pointi] > 1 - small)
                        {
                            points0[pointi] =
                                transform0(bodyMeshes_[0].bodyID_).inv()
                               .transformPoint(points[pointi]);
                        }
                        // Slerp septernion interpolation
                        else
                        {
                            points0[pointi] =
                                slerp
                                (
                                    septernion::I,
                                    septernion
                                    (
                                        transform0(bodyMeshes_[0].bodyID_)
                                    ),
                                    bodyMeshes_[0].weight_[pointi]
                                ).invTransformPoint(points[pointi]);
                        }
                    }
                    else
                    {
                        const List<septernion> transforms0(this->transforms0());
                        List<scalar> w(transforms0.size());

                        forAll(points0, pointi)
                        {
                            points0[pointi] =
                                average(transforms0, weights(pointi, w))
                               .invTransformPoint(points0[pointi]);
                        }
                    }
                }
            }
        }
    }

    points0_.transfer(points0);

    // points0 changed - set to write and check-in to database
    points0_.rename("points0");
    points0_.writeOpt() = IOobject::AUTO_WRITE;
    points0_.instance() = mesh().time().name();
    points0_.checkIn();
}


bool Foam::rigidBodyMeshMotion::write() const
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
        dict.regIOobject::writeObject
        (
            IOstream::ASCII,
            IOstream::currentVersion,
            mesh().time().writeCompression(),
            true
        )
     && displacementMotionSolver::write();
}


// ************************************************************************* //
