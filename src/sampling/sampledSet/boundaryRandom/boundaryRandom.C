/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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

#include "boundaryRandom.H"
#include "sampledSet.H"
#include "meshSearch.H"
#include "DynamicList.H"
#include "polyMesh.H"
#include "addToRunTimeSelectionTable.H"
#include "word.H"
#include "Random.H"
#include "SubField.H"
#include "barycentric2D.H"
#include "triPointRef.H"
#include "tetIndices.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace sampledSets
{
    defineTypeNameAndDebug(boundaryRandom, 0);
    addToRunTimeSelectionTable(sampledSet, boundaryRandom, word);
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::sampledSets::boundaryRandom::calcSamples
(
    DynamicList<point>& samplingPts,
    DynamicList<label>& samplingCells,
    DynamicList<label>& samplingFaces,
    DynamicList<label>& samplingSegments,
    DynamicList<scalar>& samplingCurveDist
) const
{
    // Get the patch IDs
    const labelList patchIDs(mesh().boundaryMesh().patchSet(patches_).toc());

    // Triangulate the patch faces
    DynamicList<label> triFaces, triTetPts;
    forAll(patchIDs, patchi)
    {
        const polyPatch& patch = mesh().boundaryMesh()[patchIDs[patchi]];

        forAll(patch, patchFacei)
        {
            const face& f = patch[patchFacei];
            const label facei = patchFacei + patch.start();

            for (label tetPti = 1; tetPti < f.size() - 1; ++ tetPti)
            {
                triFaces.append(facei);
                triTetPts.append(tetPti);
            }
        }
    }

    // Generate the fractions which select the processor, patch and triangle
    scalarField trisFraction(triFaces.size() + 1, 0);
    forAll(triFaces, trii)
    {
        const tetIndices tetIs
        (
            mesh().faceOwner()[triFaces[trii]],
            triFaces[trii],
            triTetPts[trii]
        );

        trisFraction[trii + 1] =
            trisFraction[trii] + tetIs.faceTri(mesh()).mag();
    }

    scalarField procsFraction(Pstream::nProcs() + 1, 0);
    {
        scalarField procsArea(Pstream::nProcs(), 0);
        procsArea[Pstream::myProcNo()] = trisFraction.last();
        Pstream::listCombineGather(procsArea, maxEqOp<scalar>());
        Pstream::listCombineScatter(procsArea);
        for(label proci = 0; proci < Pstream::nProcs(); ++ proci)
        {
            procsFraction[proci + 1] = procsFraction[proci] + procsArea[proci];
        }
    }

    if (triFaces.size())
    {
        trisFraction /= trisFraction.last();
    }

    if (procsFraction.last() != 0)
    {
        procsFraction /= procsFraction.last();
    }

    // Generate the samples
    Random rndGen(261782);
    const label proci = Pstream::myProcNo();
    for (label i = 0; i < nPoints_; ++ i)
    {
        // Request all random numbers simultaneously on all processors so that
        // the generator state stays consistent

        const scalar rProc = rndGen.scalar01();
        const scalar rTri = rndGen.scalar01();
        const barycentric2D r2D = barycentric2D01(rndGen);

        if (procsFraction[proci] < rProc && rProc <= procsFraction[proci + 1])
        {
            label trii = 0;
            while (rTri > trisFraction[trii + 1])
            {
                ++ trii;
            }

            const tetIndices tetIs
            (
                mesh().faceOwner()[triFaces[trii]],
                triFaces[trii],
                triTetPts[trii]
            );

            const barycentric r3D
            (
                rootSmall,
                (1 - rootSmall)*r2D.a(),
                (1 - rootSmall)*r2D.b(),
                (1 - rootSmall)*r2D.c()
            );

            samplingPts.append(tetIs.tet(mesh()).barycentricToPoint(r3D));
            samplingCells.append(tetIs.cell());
            samplingFaces.append(tetIs.face());
            samplingSegments.append(0);
            samplingCurveDist.append(scalar(i));
        }
    }
}


void Foam::sampledSets::boundaryRandom::genSamples()
{
    // Storage for sample points
    DynamicList<point> samplingPts;
    DynamicList<label> samplingCells;
    DynamicList<label> samplingFaces;
    DynamicList<label> samplingSegments;
    DynamicList<scalar> samplingCurveDist;

    calcSamples
    (
        samplingPts,
        samplingCells,
        samplingFaces,
        samplingSegments,
        samplingCurveDist
    );

    samplingPts.shrink();
    samplingCells.shrink();
    samplingFaces.shrink();
    samplingSegments.shrink();
    samplingCurveDist.shrink();

    setSamples
    (
        samplingPts,
        samplingCells,
        samplingFaces,
        samplingSegments,
        samplingCurveDist
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::sampledSets::boundaryRandom::boundaryRandom
(
    const word& name,
    const polyMesh& mesh,
    const meshSearch& searchEngine,
    const dictionary& dict
)
:
    sampledSet(name, mesh, searchEngine, dict),
    patches_(dict.lookup("patches")),
    nPoints_(dict.lookup<label>("nPoints"))
{
    genSamples();

    if (debug)
    {
        write(Info);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::sampledSets::boundaryRandom::~boundaryRandom()
{}


// ************************************************************************* //
