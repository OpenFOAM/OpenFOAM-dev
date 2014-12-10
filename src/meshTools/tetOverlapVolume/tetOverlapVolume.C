/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012-2013 OpenFOAM Foundation
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

#include "tetOverlapVolume.H"
#include "tetrahedron.H"
#include "tetPoints.H"
#include "polyMesh.H"
#include "OFstream.H"
#include "treeBoundBox.H"
#include "indexedOctree.H"
#include "treeDataCell.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
defineTypeNameAndDebug(tetOverlapVolume, 0);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::tetOverlapVolume::tetOverlapVolume()
{}


// * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::tetOverlapVolume::tetTetOverlapVol
(
    const tetPoints& tetA,
    const tetPoints& tetB
) const
{
    static tetPointRef::tetIntersectionList insideTets;
    label nInside = 0;
    static tetPointRef::tetIntersectionList cutInsideTets;
    label nCutInside = 0;

    tetPointRef::storeOp inside(insideTets, nInside);
    tetPointRef::storeOp cutInside(cutInsideTets, nCutInside);
    tetPointRef::sumVolOp volInside;
    tetPointRef::dummyOp outside;

    if ((tetA.tet().mag() < SMALL) || (tetB.tet().mag() < SMALL))
    {
        return 0.0;
    }

    // face0
    plane pl0(tetB[1], tetB[3], tetB[2]);
    tetA.tet().sliceWithPlane(pl0, cutInside, outside);
    if (nCutInside == 0)
    {
        return 0.0;
    }

    // face1
    plane pl1(tetB[0], tetB[2], tetB[3]);
    nInside = 0;
    for (label i = 0; i < nCutInside; i++)
    {
        const tetPointRef t = cutInsideTets[i].tet();
        t.sliceWithPlane(pl1, inside, outside);
    }
    if (nInside == 0)
    {
        return 0.0;
    }

    // face2
    plane pl2(tetB[0], tetB[3], tetB[1]);
    nCutInside = 0;
    for (label i = 0; i < nInside; i++)
    {
        const tetPointRef t = insideTets[i].tet();
        t.sliceWithPlane(pl2, cutInside, outside);
    }
    if (nCutInside == 0)
    {
        return 0.0;
    }

    // face3
    plane pl3(tetB[0], tetB[1], tetB[2]);
    for (label i = 0; i < nCutInside; i++)
    {
        const tetPointRef t = cutInsideTets[i].tet();
        t.sliceWithPlane(pl3, volInside, outside);
    }

    return volInside.vol_;
}


Foam::treeBoundBox Foam::tetOverlapVolume::pyrBb
(
    const pointField& points,
    const face& f,
    const point& fc
) const
{
    treeBoundBox bb(fc, fc);
    forAll(f, fp)
    {
        const point& pt = points[f[fp]];
        bb.min() = min(bb.min(), pt);
        bb.max() = max(bb.max(), pt);
    }
    return bb;
}


// * * * * * * * * * * * Public Member Functions  * * * * * * * * * * * * * //

bool Foam::tetOverlapVolume::cellCellOverlapMinDecomp
(
    const primitiveMesh& meshA,
    const label cellAI,
    const primitiveMesh& meshB,
    const label cellBI,
    const treeBoundBox& cellBbB,
    const scalar threshold
) const
{
    const cell& cFacesA = meshA.cells()[cellAI];
    const point& ccA = meshA.cellCentres()[cellAI];

    const cell& cFacesB = meshB.cells()[cellBI];
    const point& ccB = meshB.cellCentres()[cellBI];

    scalar vol = 0.0;

    forAll(cFacesA, cFA)
    {
        label faceAI = cFacesA[cFA];

        const face& fA = meshA.faces()[faceAI];
        const treeBoundBox pyrA = pyrBb(meshA.points(), fA, ccA);
        if (!pyrA.overlaps(cellBbB))
        {
            continue;
        }

        bool ownA = (meshA.faceOwner()[faceAI] == cellAI);

        label tetBasePtAI = 0;

        const point& tetBasePtA = meshA.points()[fA[tetBasePtAI]];

        for (label tetPtI = 1; tetPtI < fA.size() - 1; tetPtI++)
        {
            label facePtAI = (tetPtI + tetBasePtAI) % fA.size();
            label otherFacePtAI = fA.fcIndex(facePtAI);

            label pt0I = -1;
            label pt1I = -1;

            if (ownA)
            {
                pt0I = fA[facePtAI];
                pt1I = fA[otherFacePtAI];
            }
            else
            {
                pt0I = fA[otherFacePtAI];
                pt1I = fA[facePtAI];
            }

            const tetPoints tetA
            (
                ccA,
                tetBasePtA,
                meshA.points()[pt0I],
                meshA.points()[pt1I]
            );
            const treeBoundBox tetABb(tetA.bounds());


            // Loop over tets of cellB
            forAll(cFacesB, cFB)
            {
                label faceBI = cFacesB[cFB];

                const face& fB = meshB.faces()[faceBI];
                const treeBoundBox pyrB = pyrBb(meshB.points(), fB, ccB);
                if (!pyrB.overlaps(pyrA))
                {
                    continue;
                }

                bool ownB = (meshB.faceOwner()[faceBI] == cellBI);

                label tetBasePtBI = 0;

                const point& tetBasePtB = meshB.points()[fB[tetBasePtBI]];

                for (label tetPtI = 1; tetPtI < fB.size() - 1; tetPtI++)
                {
                    label facePtBI = (tetPtI + tetBasePtBI) % fB.size();
                    label otherFacePtBI = fB.fcIndex(facePtBI);

                    label pt0I = -1;
                    label pt1I = -1;

                    if (ownB)
                    {
                        pt0I = fB[facePtBI];
                        pt1I = fB[otherFacePtBI];
                    }
                    else
                    {
                        pt0I = fB[otherFacePtBI];
                        pt1I = fB[facePtBI];
                    }

                    const tetPoints tetB
                    (
                        ccB,
                        tetBasePtB,
                        meshB.points()[pt0I],
                        meshB.points()[pt1I]
                    );

                    if (!tetB.bounds().overlaps(tetABb))
                    {
                        continue;
                    }

                    vol += tetTetOverlapVol(tetA, tetB);

                    if (vol > threshold)
                    {
                        return true;
                    }
                }
            }
        }
    }

    return false;
}


Foam::scalar Foam::tetOverlapVolume::cellCellOverlapVolumeMinDecomp
(
    const primitiveMesh& meshA,
    const label cellAI,

    const primitiveMesh& meshB,
    const label cellBI,
    const treeBoundBox& cellBbB
) const
{
    const cell& cFacesA = meshA.cells()[cellAI];
    const point& ccA = meshA.cellCentres()[cellAI];

    const cell& cFacesB = meshB.cells()[cellBI];
    const point& ccB = meshB.cellCentres()[cellBI];

    scalar vol = 0.0;

    forAll(cFacesA, cFA)
    {
        label faceAI = cFacesA[cFA];

        const face& fA = meshA.faces()[faceAI];
        const treeBoundBox pyrA = pyrBb(meshA.points(), fA, ccA);
        if (!pyrA.overlaps(cellBbB))
        {
            continue;
        }

        bool ownA = (meshA.faceOwner()[faceAI] == cellAI);

        label tetBasePtAI = 0;

        const point& tetBasePtA = meshA.points()[fA[tetBasePtAI]];

        for (label tetPtI = 1; tetPtI < fA.size() - 1; tetPtI++)
        {
            label facePtAI = (tetPtI + tetBasePtAI) % fA.size();
            label otherFacePtAI = fA.fcIndex(facePtAI);

            label pt0I = -1;
            label pt1I = -1;

            if (ownA)
            {
                pt0I = fA[facePtAI];
                pt1I = fA[otherFacePtAI];
            }
            else
            {
                pt0I = fA[otherFacePtAI];
                pt1I = fA[facePtAI];
            }

            const tetPoints tetA
            (
                ccA,
                tetBasePtA,
                meshA.points()[pt0I],
                meshA.points()[pt1I]
            );
            const treeBoundBox tetABb(tetA.bounds());


            // Loop over tets of cellB
            forAll(cFacesB, cFB)
            {
                label faceBI = cFacesB[cFB];

                const face& fB = meshB.faces()[faceBI];
                const treeBoundBox pyrB = pyrBb(meshB.points(), fB, ccB);
                if (!pyrB.overlaps(pyrA))
                {
                    continue;
                }

                bool ownB = (meshB.faceOwner()[faceBI] == cellBI);

                label tetBasePtBI = 0;

                const point& tetBasePtB = meshB.points()[fB[tetBasePtBI]];

                for (label tetPtI = 1; tetPtI < fB.size() - 1; tetPtI++)
                {
                    label facePtBI = (tetPtI + tetBasePtBI) % fB.size();
                    label otherFacePtBI = fB.fcIndex(facePtBI);

                    label pt0I = -1;
                    label pt1I = -1;

                    if (ownB)
                    {
                        pt0I = fB[facePtBI];
                        pt1I = fB[otherFacePtBI];
                    }
                    else
                    {
                        pt0I = fB[otherFacePtBI];
                        pt1I = fB[facePtBI];
                    }

                    const tetPoints tetB
                    (
                        ccB,
                        tetBasePtB,
                        meshB.points()[pt0I],
                        meshB.points()[pt1I]
                    );
                    if (!tetB.bounds().overlaps(tetABb))
                    {
                        continue;
                    }

                    vol += tetTetOverlapVol(tetA, tetB);
                }
            }
        }
    }

    return vol;
}


Foam::labelList Foam::tetOverlapVolume::overlappingCells
(
    const polyMesh& fromMesh,
    const polyMesh& toMesh,
    const label iTo
) const
{
    const indexedOctree<treeDataCell>& treeA = fromMesh.cellTree();

    treeBoundBox bbB(toMesh.points(), toMesh.cellPoints()[iTo]);

    return treeA.findBox(bbB);
}


// ************************************************************************* //
