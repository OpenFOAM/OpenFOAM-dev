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

#include "multiRigidBody_pointMeshMover.H"
#include "polyTopoChangeMap.H"
#include "pointDist.H"
#include "pointConstraints.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace pointMeshMovers
{
    defineTypeNameAndDebug(multiRigidBody, 0);
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::List<Foam::scalar>& Foam::pointMeshMovers::multiRigidBody::weights
(
    const label pointi,
    List<scalar>& w
) const
{
    // Initialise to 1 for the far-field weight
    scalar sumw = 1;

    // Accumulate the weighted body and exterior weights
    forAll(bodyMeshes_, bi)
    {
        const scalar wbi = bodyMeshes_[bi].weight_[pointi];
        w[bi] = wbi/pow(max(1 - wbi, small), 0.62);
        sumw += w[bi];
    }

    // Calculate the limiter for wbi/(1 - wbi) to ensure the sum(wbi) = 1
    const scalar lambda = 1/sumw;

    // Sum the limited body weights except the exterior weight
    sumw = 0;
    for(label bi=0; bi<bodyMeshes_.size() - 1; bi++)
    {
        w[bi] = lambda*w[bi];
        sumw += w[bi];
    }

    // Calculate the exterior weight
    w[bodyMeshes_.size() - 1] = 1 - sumw;

    return w;
}


void Foam::pointMeshMovers::multiRigidBody::bodyMesh::calcWeights
(
    const pointField& points0
)
{
    if (do_ > 0)
    {
        const pointMesh& pMesh = pointMesh::New(mesh_);

        // Calculate scaling factor everywhere for body
        const pointDist pDist
        (
            pMesh,
            patchSet_,
            pointZoneSet_,
            labelHashSet(),
            labelHashSet(),
            points0,
            do_
        );

        weight_.primitiveFieldRef() = weight(pDist.primitiveField());

        pointConstraints::New(pMesh).constrain(weight_);
    }
    else
    {
        weight_.primitiveFieldRef() = 0;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::pointMeshMovers::multiRigidBody::bodyMesh::bodyMesh
(
    const polyMesh& mesh,
    const word& name,
    const dictionary& dict
)
:
    mesh_(mesh),
    name_(name),
    patches_(dict.lookupOrDefault("patches", wordReList::null())),
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


Foam::pointMeshMovers::multiRigidBody::bodyMesh::bodyMesh
(
    const polyMesh& mesh,
    const word& name
)
:
    mesh_(mesh),
    name_(name),
    di_(0),
    do_(0),
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


Foam::pointMeshMovers::multiRigidBody::multiRigidBody
(
    const polyMesh& mesh,
    const dictionary& dict
)
:
    displacementPoints0(mesh, dict, typeName)
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
        bodyMeshes_.append(new bodyMesh(mesh, "body", dict));
    }

    // Append the body corresponding to the outer boundary of the region
    if (dict.isDict("exterior"))
    {
        bodyMeshes_.append
        (
            new bodyMesh
            (
                mesh,
                "exterior",
                dict.subDict("exterior")
            )
        );
    }
    else
    {
        bodyMeshes_.append
        (
            new bodyMesh(mesh, "exterior")
        );
    }

    // If the patches of the external body are not specified
    // include all non-constraint patches except those of the moving bodies
    if (!bodyMeshes_.last().patchSet_.size())
    {
        labelHashSet& externalPatches = bodyMeshes_.last().patchSet_;

        forAll(mesh.boundary(), patchi)
        {
            if (!mesh.boundary()[patchi].constraint())
            {
                externalPatches.insert(patchi);
            }
        }

        for(label bi = 0; bi < bodyMeshes_.size() - 1; bi++)
        {
            externalPatches -= bodyMeshes_[bi].patchSet_;
        }
    }

    // Calculate scaling factor everywhere for each meshed body
    // from the current mesh points
    forAll(bodyMeshes_, bi)
    {
        // bodyMeshes_[bi].calcWeights(points0());
        bodyMeshes_[bi].calcWeights(mesh.points());
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::pointMeshMovers::multiRigidBody::~multiRigidBody()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Type Foam::pointMeshMovers::multiRigidBody::bodyMesh::weight
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
Foam::pointMeshMovers::multiRigidBody::newPoints()
{
    if (poly().nPoints() != points0().size())
    {
        FatalErrorInFunction
            << "The number of points in the mesh seems to have changed." << endl
            << "In constant/polyMesh there are " << points0().size()
            << " points; in the current mesh there are " << poly().nPoints()
            << " points." << exit(FatalError);
    }

    moveBodies();

    const pointField& points0 = this->points0();
    tmp<pointField> tpoints(new pointField(points0.size()));
    pointField& points(tpoints.ref());

    // Calculate the transformations of all the bodies and the exterior
    const List<septernion> transforms0(this->transforms0());

    // Storage for the body weights
    List<scalar> w(transforms0.size());

    // Transform the points with the SLERP average of the body transformations
    forAll(points0, pointi)
    {
        points[pointi] =
            average(transforms0, weights(pointi, w))
           .transformPoint(points0[pointi]);
    }

    return tpoints;
}


void Foam::pointMeshMovers::multiRigidBody::topoChange
(
    const polyTopoChangeMap& map
)
{
    // Get the new points either from the map or the mesh
    const pointField& points = poly().points();

    const pointMesh& pMesh = pointMesh::New(poly());

    pointField newPoints0(poly().points());

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
            if (bodyMeshes_[bi].do_ > 0)
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

                scalarField& weight = bodyMeshes_[bi].weight_;

                forAll(newPoints0, pointi)
                {
                    const label oldPointi = map.pointMap()[pointi];

                    if (map.reversePointMap()[oldPointi] != pointi)
                    {
                        weight[pointi] = bodyMeshes_[bi].weight(pDist[pointi]);
                    }
                }

                pointConstraints::New(pMesh).constrain(bodyMeshes_[bi].weight_);
            }
            else
            {
                bodyMeshes_[bi].weight_ = 0;
            }
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

        // Calculate the transformations of all the bodies and the exterior
        const List<septernion> transforms0(this->transforms0());

        // Storage for the body weights
        List<scalar> w(transforms0.size());

        // Inverse Transform the points with the SLERP average
        // of the body transformations
        forAll(newPoints0, pointi)
        {
            const label oldPointi = map.pointMap()[pointi];

            if (map.reversePointMap()[oldPointi] != pointi)
            {
                newPoints0[pointi] =
                    average(transforms0, weights(pointi, w))
                   .invTransformPoint(points[pointi]);
            }
        }
    }

    // Move into base class storage
    points0_.primitiveFieldRef() = newPoints0;

    // Mark the changed points0 to be written automatically
    points0_.writeOpt() = IOobject::AUTO_WRITE;
    points0_.instance() = poly().time().name();
}


void Foam::pointMeshMovers::multiRigidBody::mapMesh(const polyMeshMap& map)
{
    displacementPoints0::mapMesh(map);

    pointField& points0 = this->points0();

    // Calculate scaling factor everywhere for each meshed body
    forAll(bodyMeshes_, bi)
    {
        bodyMeshes_[bi].calcWeights(points0);
    }

    // Calculate the transformations of all the bodies and the exterior
    const List<septernion> transforms0(this->transforms0());

    // Storage for the body weights
    List<scalar> w(transforms0.size());

    // Inverse Transform the points0 with the SLERP average
    // of the body transformations
    forAll(points0, pointi)
    {
        points0[pointi] =
            average(transforms0, weights(pointi, w))
           .invTransformPoint(points0[pointi]);
    }

    // Mark the changed points0 to be written automatically
    points0_.writeOpt() = IOobject::AUTO_WRITE;
    points0_.instance() = poly().time().name();
}


// ************************************************************************* //
