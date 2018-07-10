/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2015-2018 OpenFOAM Foundation
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

#include "blockMesh.H"

// * * * * * * * * * * * * * * * Private Functions * * * * * * * * * * * * * //

namespace Foam
{

// faces
// 6
// (
//     4(0 4 7 3)   // x-min
//     4(1 2 6 5)   // x-max
//     4(0 1 5 4)   // y-min
//     4(3 7 6 2)   // y-max
//     4(0 3 2 1)   // z-min
//     4(4 5 6 7)   // z-max
// );

// Face-edge directions
static const int faceEdgeDirs[6][4] =
{
    {2, 1, -2, -1},
    {1, 2, -1, -2},
    {1, 2, -1, -2},
    {2, 1, -2, -1},
    {2, 1, -2, -1},
    {1, 2, -1, -2}
};

// The face-face-rotation direction correspondence map
static Pair<int> faceFaceRotMap[6][6][4];

// Generate the face-face-rotation direction correspondence map
void genFaceFaceRotMap()
{
    for(int facePi=0; facePi<6; facePi++)
    {
        for(int faceNi=0; faceNi<6; faceNi++)
        {
            for(int rot=0; rot<4; rot++)
            {
                Pair<int>& map = faceFaceRotMap[facePi][faceNi][rot];

                for(int Pp=0; Pp<2; Pp++)
                {
                    int Pdir = faceEdgeDirs[facePi][Pp];
                    int Np = (3 - Pp + rot)%4;
                    int Ndir = faceEdgeDirs[faceNi][Np];
                    map[Pdir-1] = -Ndir;
                }

                // Handle sign change due to match-face transpose
                if (mag(map[0]) == 2 && map[0]*map[1] < 0)
                {
                    map[0] = -map[0];
                    map[1] = -map[1];
                }
            }
        }
    }
}

// Return the direction map for the merge-faces
Pair<int> faceMap
(
    const label facePi,
    const face& faceP,
    const label faceNi,
    const face& faceN
)
{
    // Search for the point on faceN corresponding to the 0-point on faceP
    for(int rot=0; rot<4; rot++)
    {
        if (faceN[rot] == faceP[0])
        {
            return faceFaceRotMap[facePi][faceNi][rot];
        }
    }

    FatalErrorInFunction
        << "Cannot find point correspondence for faces "
        << faceP << " and " << faceN
        << exit(FatalError);

    return Pair<int>(0, 0);
}

// Set the block and face indices for all the merge faces
void setBlockFaceCorrespondence
(
    const cellList& topoCells,
    const faceList::subList& topoInternalFaces,
    const labelList& topoFaceCell,
    List<Pair<label>>& mergeBlock
)
{
    forAll(topoInternalFaces, topoFacei)
    {
        label topoPi = topoFaceCell[topoFacei];
        const labelList& topoPfaces = topoCells[topoPi];

        bool foundFace = false;
        label topoPfacei;
        for
        (
            topoPfacei = 0;
            topoPfacei < topoPfaces.size();
            topoPfacei++
        )
        {
            if (topoPfaces[topoPfacei] == topoFacei)
            {
                foundFace = true;
                break;
            }
        }

        if (!foundFace)
        {
            FatalErrorInFunction
                << "Cannot find merge face for block " << topoPi
                << exit(FatalError);
        }

        mergeBlock[topoFacei].first() = topoPi;
        mergeBlock[topoFacei].second() = topoPfacei;
    }
}

// Return the number of divisions in each direction for the face
Pair<label> faceNij(const label facei, const block& block)
{
    Pair<label> fnij;

    int i = facei/2;

    if (i == 0)
    {
        fnij.first() = block.density().y() + 1;
        fnij.second() = block.density().z() + 1;
    }
    else if (i == 1)
    {
        fnij.first() = block.density().x() + 1;
        fnij.second() = block.density().z() + 1;
    }
    else if (i == 2)
    {
        fnij.first() = block.density().x() + 1;
        fnij.second() = block.density().y() + 1;
    }

    return fnij;
}

// Sign the index corresponding to the map
inline label signIndex(const int map, const label i)
{
    return map < 0 ? -i-1 : i;
}

// Reverse a signed index with the number of divisions
inline label unsignIndex(const label i, const label ni)
{
    return i >= 0 ? i : ni + i + 1;
}

// Return the mapped index
inline label mapij(const int map, const label i, const label j)
{
    return signIndex(map, mag(map) == 1 ? i : j);
}

// Return the face point index
inline label facePoint
(
    const int facei,
    const block& block,
    const label i,
    const label j
)
{
    switch (facei)
    {
        case 0:
            return block.pointLabel(0, i, j);
        case 1:
            return block.pointLabel(block.density().x(), i, j);
        case 2:
            return block.pointLabel(i, 0, j);
        case 3:
            return block.pointLabel(i, block.density().y(), j);
        case 4:
            return block.pointLabel(i, j, 0);
        case 5:
            return block.pointLabel(i, j, block.density().z());
        default:
            return -1;
    }
}

// Return the neighbour face point from the signed indices
inline label facePointN
(
    const block& block,
    const label i,
    const label j,
    const label k
)
{
    return block.pointLabel
    (
        unsignIndex(i, block.density().x()),
        unsignIndex(j, block.density().y()),
        unsignIndex(k, block.density().z())
    );
}

// Return the neighbour face point from the mapped indices
inline label facePointN
(
    const int facei,
    const block& block,
    const label i,
    const label j
)
{
    switch (facei)
    {
        case 0:
            return facePointN(block, 0, i, j);
        case 1:
            return facePointN(block, block.density().x(), i, j);
        case 2:
            return facePointN(block, i, 0, j);
        case 3:
            return facePointN(block, i, block.density().y(), j);
        case 4:
            return facePointN(block, i, j, 0);
        case 5:
            return facePointN(block, i, j, block.density().z());
        default:
            return -1;
    }
}

// Return the neighbour face point using the map
inline label facePointN
(
    const int facei,
    const Pair<int>& fmap,
    const block& block,
    const label i,
    const label j
)
{
    return facePointN(facei, block, mapij(fmap[0], i, j), mapij(fmap[1], i, j));
}

}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::blockMesh::calcMergeInfoFast()
{
    // Generate the static face-face map
    genFaceFaceRotMap();

    const blockList& blocks = *this;

    if (verboseOutput)
    {
        Info<< "Creating block offsets" << endl;
    }

    blockOffsets_.setSize(blocks.size());

    nPoints_ = 0;
    nCells_  = 0;

    forAll(blocks, blocki)
    {
        blockOffsets_[blocki] = nPoints_;

        nPoints_ += blocks[blocki].nPoints();
        nCells_  += blocks[blocki].nCells();
    }

    if (verboseOutput)
    {
        Info<< "Creating merge list using the fast topological search"
            << flush;
    }

    // Size merge list and initialize to -1
    mergeList_.setSize(nPoints_, -1);

    // Block mesh topology
    const pointField& topoPoints = topology().points();
    const cellList& topoCells = topology().cells();
    const faceList& topoFaces = topology().faces();
    const labelList& topoFaceOwn = topology().faceOwner();
    const labelList& topoFaceNei = topology().faceNeighbour();

    // Topological merging is only necessary for internal block faces
    // Note edge and face collapse may apply to boundary faces
    // but is not yet supported in the "fast" algorithm
    const faceList::subList topoInternalFaces
    (
        topoFaces,
        topology().nInternalFaces()
    );

    List<Pair<label>> mergeBlockP(topoInternalFaces.size());
    setBlockFaceCorrespondence
    (
        topoCells,
        topoInternalFaces,
        topoFaceOwn,
        mergeBlockP
    );

    List<Pair<label>> mergeBlockN(topoInternalFaces.size());
    setBlockFaceCorrespondence
    (
        topoCells,
        topoInternalFaces,
        topoFaceNei,
        mergeBlockN
    );

    if (debug)
    {
        Info<< endl;
    }

    forAll(topoInternalFaces, topoFacei)
    {
        if (debug)
        {
            Info<< "Processing face " << topoFacei << endl;
        }

        label blockPi = mergeBlockP[topoFacei].first();
        label blockPfacei = mergeBlockP[topoFacei].second();

        label blockNi = mergeBlockN[topoFacei].first();
        label blockNfacei = mergeBlockN[topoFacei].second();

        Pair<int> fmap
        (
            faceMap
            (
                blockPfacei,
                blocks[blockPi].blockShape().faces()[blockPfacei],
                blockNfacei,
                blocks[blockNi].blockShape().faces()[blockNfacei]
            )
        );

        if (debug)
        {
            Info<< "    Face map for faces "
                << blocks[blockPi].blockShape().faces()[blockPfacei] << " "
                << blocks[blockNi].blockShape().faces()[blockNfacei] << ": "
                << fmap << endl;
        }

        const pointField& blockPpoints = blocks[blockPi].points();
        const pointField& blockNpoints = blocks[blockNi].points();

        Pair<label> Pnij(faceNij(blockPfacei, blocks[blockPi]));

        // Check block subdivision correspondence
        {
            Pair<label> Nnij(faceNij(blockNfacei, blocks[blockNi]));
            Pair<label> NPnij;
            NPnij[0] = Nnij[mag(fmap[0]) - 1];
            NPnij[1] = Nnij[mag(fmap[1]) - 1];

            if (Pnij != NPnij)
            {
                FatalErrorInFunction
                    << "Sub-division mismatch between face "
                    << blockPfacei << " of block " << blockPi << Pnij
                    << " and face "
                    << blockNfacei << " of block " << blockNi << Nnij
                    << exit(FatalError);
            }
        }

        // Calculate a suitable test distance from the bounding box of the face.
        // Note this is used only as a sanity check and for diagnostics in
        // case there is a grading inconsistency.
        const boundBox bb(topoCells[blockPi].points(topoFaces, topoPoints));
        const scalar testSqrDist = magSqr(1e-6*bb.span());

        // Accumulate the maximum merge distance for diagnostics
        scalar maxSqrDist = 0;

        for (label j=0; j<Pnij.second(); j++)
        {
            for (label i=0; i<Pnij.first(); i++)
            {
                label blockPpointi =
                    facePoint(blockPfacei, blocks[blockPi], i, j);

                label blockNpointi =
                    facePointN(blockNfacei, fmap, blocks[blockNi], i, j);

                scalar sqrDist
                (
                    magSqr
                    (
                        blockPpoints[blockPpointi]
                      - blockNpoints[blockNpointi]
                    )
                );

                if (sqrDist > testSqrDist)
                {
                    FatalErrorInFunction
                        << "Point merge failure between face "
                        << blockPfacei << " of block " << blockPi
                        << " and face "
                        << blockNfacei << " of block " << blockNi
                        << endl
                        << "    Points: " << blockPpoints[blockPpointi]
                        << " " << blockNpoints[blockNpointi]
                        << endl
                        << "    This may be due to inconsistent grading."
                        << exit(FatalError);
                }

                maxSqrDist = max(maxSqrDist, sqrDist);

                label Ppointi = blockPpointi + blockOffsets_[blockPi];
                label Npointi = blockNpointi + blockOffsets_[blockNi];

                label minPNi = min(Ppointi, Npointi);

                if (mergeList_[Ppointi] != -1)
                {
                    minPNi = min(minPNi, mergeList_[Ppointi]);
                }

                if (mergeList_[Npointi] != -1)
                {
                    minPNi = min(minPNi, mergeList_[Npointi]);
                }

                mergeList_[Ppointi] = mergeList_[Npointi] = minPNi;
            }
        }

        if (debug)
        {
            Info<< "    Max distance between merge points: "
                << sqrt(maxSqrDist) << endl;
        }
    }


    bool changedPointMerge = false;
    label nPasses = 0;

    do
    {
        changedPointMerge = false;
        nPasses++;

        forAll(topoInternalFaces, topoFacei)
        {
            label blockPi = mergeBlockP[topoFacei].first();
            label blockPfacei = mergeBlockP[topoFacei].second();

            label blockNi = mergeBlockN[topoFacei].first();
            label blockNfacei = mergeBlockN[topoFacei].second();

            Pair<int> fmap
            (
                faceMap
                (
                    blockPfacei,
                    blocks[blockPi].blockShape().faces()[blockPfacei],
                    blockNfacei,
                    blocks[blockNi].blockShape().faces()[blockNfacei]
                )
            );

            Pair<label> Pnij(faceNij(blockPfacei, blocks[blockPi]));

            for (label j=0; j<Pnij.second(); j++)
            {
                for (label i=0; i<Pnij.first(); i++)
                {
                    label blockPpointi =
                        facePoint(blockPfacei, blocks[blockPi], i, j);

                    label blockNpointi =
                        facePointN(blockNfacei, fmap, blocks[blockNi], i, j);

                    label Ppointi =
                        blockPpointi + blockOffsets_[blockPi];

                    label Npointi =
                        blockNpointi + blockOffsets_[blockNi];

                    if (mergeList_[Ppointi] != mergeList_[Npointi])
                    {
                        changedPointMerge = true;

                        mergeList_[Ppointi]
                          = mergeList_[Npointi]
                          = min(mergeList_[Ppointi], mergeList_[Npointi]);
                    }
                }
            }
        }

        if (verboseOutput)
        {
            Info<< "." << flush;
        }

        if (nPasses > 100)
        {
            FatalErrorInFunction
                << "Point merging failed after 100 passes."
                << exit(FatalError);
        }

    } while (changedPointMerge);


    // Sort merge list and count number of unique points
    label nUniqPoints = 0;

    forAll(mergeList_, pointi)
    {
        if (mergeList_[pointi] > pointi)
        {
            FatalErrorInFunction
                << "Merge list contains point index out of range"
                << exit(FatalError);
        }

        if (mergeList_[pointi] == -1 || mergeList_[pointi] == pointi)
        {
            mergeList_[pointi] = nUniqPoints++;
        }
        else
        {
            mergeList_[pointi] = mergeList_[mergeList_[pointi]];
        }
    }

    // Correct number of points in mesh
    nPoints_ = nUniqPoints;
}


// ************************************************************************* //
