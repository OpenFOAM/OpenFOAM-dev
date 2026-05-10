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

#include "multiRigidBodyMeshMotion.H"
#include "polyTopoChangeMap.H"
#include "pointDist.H"
#include "pointConstraints.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(multiRigidBodyMeshMotion, 0);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::List<Foam::scalar>& Foam::multiRigidBodyMeshMotion::weights
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
    const scalar lambda = 1/sum1mw;

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

Foam::multiRigidBodyMeshMotion::bodyMesh::bodyMesh
(
    const polyMesh& mesh,
    const word& name,
    const dictionary& dict
)
:
    name_(name),
    patches_(wordReList(dict.lookup("patches"))),
    patchSet_(mesh.boundary().patchSet(patches_)),
    pointZones_(dict.lookupOrDefault("pointZones", wordReList::null())),
    pointZoneSet_(mesh.pointZones().zoneSet(pointZones_)),
    di_(dict.lookupOrDefault<scalar>("innerDistance", 0.0)),
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
    ),
    bodyIndex(-1)
{}


Foam::multiRigidBodyMeshMotion::multiRigidBodyMeshMotion
(
    const word& name,
    const polyMesh& mesh,
    const dictionary& dict
)
:
    points0MotionSolver(name, mesh, dict, typeName)
{
    if (dict.isDict("bodies"))
    {
        const dictionary& bodiesDict = dict.subDict("bodies");

        forAllConstIter(IDLList<entry>, bodiesDict, iter)
        {
            const dictionary& bodyDict = iter().dict();

            if (bodyDict.found("patches"))
            {
                bodyMeshes_.append
                (
                    new bodyMesh
                    (
                        mesh,
                        iter().keyword(),
                        bodyDict
                    )
                );
            }
        }
    }
    else
    {
        bodyMeshes_.append(new bodyMesh(mesh, name, dict));
    }

    const pointMesh& pMesh = pointMesh::New(mesh);

    // Calculate scaling factor everywhere for each meshed body
    forAll(bodyMeshes_, bi)
    {
        const pointDist pDist
        (
            pMesh,
            bodyMeshes_[bi].patchSet_,
            bodyMeshes_[bi].pointZoneSet_,
            labelHashSet(),
            labelHashSet(),
            points0(),
            bodyMeshes_[bi].do_
        );

        bodyMeshes_[bi].weight_.primitiveFieldRef() =
            bodyMeshes_[bi].weight(pDist.primitiveField());

        pointConstraints::New(pMesh).constrain(bodyMeshes_[bi].weight_);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::multiRigidBodyMeshMotion::~multiRigidBodyMeshMotion()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Type Foam::multiRigidBodyMeshMotion::bodyMesh::weight
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
Foam::multiRigidBodyMeshMotion::newPoints()
{
    if (mesh().nPoints() != points0().size())
    {
        FatalErrorInFunction
            << "The number of points in the mesh seems to have changed." << endl
            << "In constant/polyMesh there are " << points0().size()
            << " points; in the current mesh there are " << mesh().nPoints()
            << " points." << exit(FatalError);
    }

    moveBodies();

    const pointField& points0 = this->points0();
    vectorField pointDisplacement(points0.size());

    // Update the displacements
    if (bodyMeshes_.size() == 1)
    {
        const septernion transform0(this->transforms0()[0]);
        const scalarField& weight = bodyMeshes_[0].weight_;

        forAll(points0, pointi)
        {
            // Don't move where weight ~= 0
            if (weight[pointi] <= small)
            {
                pointDisplacement[pointi] = Zero;
            }
            // Use solid-body motion where weight ~= 1
            else if (weight[pointi] > 1 - small)
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
    else
    {
        List<septernion> transforms0(this->transforms0());
        transforms0.append(septernion::I);
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
    // pointConstraints::New
    // (
    //     pointDisplacement.mesh()
    // ).constrainDisplacement(pointDisplacement);

    return points0 + pointDisplacement;
}


void Foam::multiRigidBodyMeshMotion::topoChange(const polyTopoChangeMap& map)
{
    // pointMesh already updates pointFields

    // Get the new points either from the map or the mesh
    const pointField& points = mesh().points();

    const pointMesh& pMesh = pointMesh::New(mesh());

    pointField newPoints0(mesh().points());

    forAll(newPoints0, pointi)
    {
        if (map.pointMap()[pointi] < 0)
        {
            FatalErrorInFunction
                << "Cannot determine co-ordinates of introduced vertices."
                << " New vertex " << pointi << " at co-ordinate "
                << points[pointi] << exit(FatalError);
        }
    }

    // Iterate to update the transformation of the new points to the
    // corresponding points0, required because the body-point weights are
    // calculated for points0
    for (int iter=0; iter<3; iter++)
    {
        // Calculate scaling factor everywhere for each meshed body
        forAll(bodyMeshes_, bi)
        {
            const pointDist pDist
            (
                pMesh,
                bodyMeshes_[bi].patchSet_,
                bodyMeshes_[bi].pointZoneSet_,
                labelHashSet(),
                labelHashSet(),
                newPoints0,
                bodyMeshes_[bi].do_
            );

            pointScalarField& weight = bodyMeshes_[bi].weight_;

            forAll(newPoints0, pointi)
            {
                const label oldPointi = map.pointMap()[pointi];

                if (map.reversePointMap()[oldPointi] != pointi)
                {
                    weight[pointi] = bodyMeshes_[bi].weight(pDist[pointi]);
                }
            }

            pointConstraints::New(pMesh).constrain(weight);
        }

        // Set directly mapped points
        forAll(newPoints0, pointi)
        {
            const label oldPointi = map.pointMap()[pointi];

            if (map.reversePointMap()[oldPointi] == pointi)
            {
                newPoints0[pointi] = points0_[oldPointi];
            }
        }

        // Interpolate indirectly mapped points
        if (bodyMeshes_.size() == 1)
        {
            const septernion transform0(this->transforms0()[0]);
            const scalarField& weight = bodyMeshes_[0].weight_;

            forAll(newPoints0, pointi)
            {
                const label oldPointi = map.pointMap()[pointi];

                if (map.reversePointMap()[oldPointi] != pointi)
                {
                    // Don't move where weight ~= 0
                    if (weight[pointi] <= small)
                    {
                        newPoints0[pointi] = points[pointi];
                    }
                    // Use solid-body motion where weight ~= 1
                    else if (weight[pointi] > 1 - small)
                    {
                        newPoints0[pointi] =
                            transform0.invTransformPoint(points[pointi]);
                    }
                    // Slerp septernion interpolation
                    else
                    {
                        newPoints0[pointi] =
                            slerp(septernion::I, transform0, weight[pointi])
                           .invTransformPoint(points[pointi]);
                    }
                }
            }
        }
        else
        {
            const List<septernion> transforms0(this->transforms0());
            List<scalar> w(transforms0.size());

            forAll(newPoints0, pointi)
            {
                const label oldPointi = map.pointMap()[pointi];

                if (map.reversePointMap()[oldPointi] != pointi)
                {
                    newPoints0[pointi] =
                        average(transforms0, weights(pointi, w))
                       .invTransformPoint(newPoints0[pointi]);
                }
            }
        }
    }

    // Move into base class storage and mark as to-be-written
    points0_.primitiveFieldRef() = newPoints0;
    points0_.writeOpt() = IOobject::AUTO_WRITE;
    points0_.instance() = mesh().time().name();
}


// ************************************************************************* //
