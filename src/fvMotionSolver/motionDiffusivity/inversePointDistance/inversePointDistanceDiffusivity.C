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

#include "inversePointDistanceDiffusivity.H"
#include "surfaceFields.H"
#include "HashSet.H"
#include "pointEdgePoint.H"
#include "PointEdgeWave.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(inversePointDistanceDiffusivity, 0);

    addToRunTimeSelectionTable
    (
        motionDiffusivity,
        inversePointDistanceDiffusivity,
        Istream
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::inversePointDistanceDiffusivity::inversePointDistanceDiffusivity
(
    const fvMesh& mesh,
    Istream& mdData
)
:
    motionDiffusivity(mesh),
    patchNames_(mdData)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::inversePointDistanceDiffusivity::~inversePointDistanceDiffusivity()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::surfaceScalarField>
Foam::inversePointDistanceDiffusivity::operator()() const
{
    tmp<surfaceScalarField> tfaceDiffusivity
    (
        new surfaceScalarField
        (
            IOobject
            (
                "faceDiffusivity",
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh(),
            dimensionedScalar(dimless, 1.0)
        )
    );

    surfaceScalarField& faceDiffusivity = tfaceDiffusivity.ref();

    const polyBoundaryMesh& bdry = mesh().boundaryMesh();

    const labelHashSet patchSet(bdry.patchSet(patchNames_));

    label nPatchEdges = 0;

    forAllConstIter(labelHashSet, patchSet, iter)
    {
        nPatchEdges += bdry[iter.key()].nEdges();
    }

    // Distance to wall on points and edges.
    List<pointEdgePoint> pointWallDist(mesh().nPoints());
    List<pointEdgePoint> edgeWallDist(mesh().nEdges());

    int dummyTrackData = 0;


    {
        // Seeds
        List<pointEdgePoint> seedInfo(nPatchEdges);
        labelList seedPoints(nPatchEdges);

        nPatchEdges = 0;

        forAllConstIter(labelHashSet, patchSet, iter)
        {
            const polyPatch& patch = bdry[iter.key()];

            const labelList& meshPoints = patch.meshPoints();

            forAll(meshPoints, i)
            {
                const label pointi = meshPoints[i];

                if (!pointWallDist[pointi].valid(dummyTrackData))
                {
                    // Not yet seeded
                    seedInfo[nPatchEdges] = pointEdgePoint
                    (
                        mesh().points()[pointi],
                        0.0
                    );
                    seedPoints[nPatchEdges] = pointi;
                    pointWallDist[pointi] = seedInfo[nPatchEdges];

                    nPatchEdges++;
                }
            }
        }
        seedInfo.setSize(nPatchEdges);
        seedPoints.setSize(nPatchEdges);

        // Do calculations
        PointEdgeWave<pointEdgePoint> waveInfo
        (
            mesh(),
            seedPoints,
            seedInfo,

            pointWallDist,
            edgeWallDist,
            mesh().globalData().nTotalPoints(),// max iterations
            dummyTrackData
        );
    }


    for (label facei=0; facei<mesh().nInternalFaces(); facei++)
    {
        const face& f = mesh().faces()[facei];

        scalar dist = 0;

        forAll(f, fp)
        {
            dist += sqrt(pointWallDist[f[fp]].distSqr());
        }
        dist /= f.size();

        faceDiffusivity[facei] = 1.0/dist;
    }


    surfaceScalarField::Boundary& faceDiffusivityBf =
        faceDiffusivity.boundaryFieldRef();

    forAll(faceDiffusivityBf, patchi)
    {
        fvsPatchScalarField& bfld = faceDiffusivityBf[patchi];

        if (patchSet.found(patchi))
        {
            const labelUList& faceCells = bfld.patch().faceCells();

            forAll(bfld, i)
            {
                const cell& ownFaces = mesh().cells()[faceCells[i]];

                labelHashSet cPoints(4*ownFaces.size());

                scalar dist = 0;

                forAll(ownFaces, ownFacei)
                {
                    const face& f = mesh().faces()[ownFaces[ownFacei]];

                    forAll(f, fp)
                    {
                        if (cPoints.insert(f[fp]))
                        {
                            dist += sqrt(pointWallDist[f[fp]].distSqr());
                        }
                    }
                }
                dist /= cPoints.size();

                bfld[i] = 1.0/dist;
            }
        }
        else
        {
            const label start = bfld.patch().start();

            forAll(bfld, i)
            {
                const face& f = mesh().faces()[start+i];

                scalar dist = 0;

                forAll(f, fp)
                {
                    dist += sqrt(pointWallDist[f[fp]].distSqr());
                }
                dist /= f.size();

                bfld[i] = 1.0/dist;
            }
        }
    }

    return tfaceDiffusivity;
}


// ************************************************************************* //
