/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2025 OpenFOAM Foundation
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

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::blockMesh::check(const polyMesh& bm, const dictionary& dict) const
{
    const pointField& points = bm.points();
    const cellList& cells = bm.cells();
    const polyPatchList& patches = bm.boundaryMesh();

    label nBoundaryFaces = 0;
    forAll(cells, celli)
    {
        nBoundaryFaces += cells[celli].nFaces();
    }

    nBoundaryFaces -= 2*bm.nInternalFaces();

    label nDefinedBoundaryFaces = 0;
    forAll(patches, patchi)
    {
        nDefinedBoundaryFaces += patches[patchi].size();
    }

    if (verboseOutput)
    {
        Info<< nl << "Basic statistics" << nl
            << "    Number of internal faces : "
            << bm.nInternalFaces() << nl
            << "    Number of boundary faces : "
            << nBoundaryFaces << nl
            << "    Number of defined boundary faces : "
            << nDefinedBoundaryFaces << nl
            << "    Number of undefined boundary faces : "
            << nBoundaryFaces - nDefinedBoundaryFaces << nl;

        if (nBoundaryFaces - nDefinedBoundaryFaces > 0)
        {
            Info<< "        (Warning : only leave undefined the front "
                << "and back planes of 2D planar geometries!)" << endl;
        }
    }

    Info<< nl << "Checking topology" << endl;

    OStringStream oss;

    // Check for duplicate curved edge definitions
    forAll(edges_, cei)
    {
        for (label cej=cei+1; cej<edges_.size(); cej++)
        {
            if (edges_[cei].compare(edges_[cej]) != 0)
            {
                oss << "    Curved edge #" << cei << ' ';
                edges_[cei].write(oss, dict);
                oss << " is a duplicate of curved edge #" << cej << ' ';
                edges_[cej].write(oss, dict);
                oss << endl;
                break;
            }
        }
    }

    // Check curved-edge/block-edge correspondence
    //
    // Loop over the edges of each block rather than the edgeList of the
    // topological mesh due to problems with calcEdges for blocks with
    // repeated point labels
    const blockList& blocks = *this;

    forAll(edges_, cei)
    {
        bool found = false;

        forAll(blocks, blocki)
        {
            edgeList edges = blocks[blocki].blockShape().edges();

            forAll(edges, ei)
            {
                found = edges_[cei].compare(edges[ei][0], edges[ei][1]) != 0;
                if (found) break;
            }
            if (found) break;
        }

        if (!found)
        {
            oss << "    Curved edge #" << cei << ' ';
            edges_[cei].write(oss, dict);
            oss << " does not correspond to a block edge" << endl;
        }
    }

    const faceList& faces = bm.faces();

    // Check for duplicate curved face definitions
    forAll(faces_, cfi)
    {
        for (label cfj=cfi+1; cfj<faces_.size(); cfj++)
        {
            if (faces_[cfi].compare(faces_[cfj]) != 0)
            {
                oss << "    Curved face #" << cfi << ' ';
                faces_[cfi].write(oss, dict);
                oss << " is a duplicate of curved face #" << cfj << ' ';
                faces_[cfj].write(oss, dict);
                oss << endl;
                break;
            }
        }
    }

    // Check curved-face/block-face correspondence
    forAll(faces_, cfi)
    {
        bool found = false;

        forAll(faces, fi)
        {
            found = faces_[cfi].compare(faces[fi]) != 0;
            if (found) break;
        }

        if (!found)
        {
            oss << "    Curved face #" << cfi << ' ';
            faces_[cfi].write(oss, dict);
            oss << " does not correspond to a block face" << endl;
        }
    }

    // Check patch-face consistency
    forAll(patches, patchi)
    {
        const faceList& Patch = patches[patchi];

        forAll(Patch, patchFacei)
        {
            const face& patchFace = Patch[patchFacei];
            bool patchFaceOK = false;

            forAll(cells, celli)
            {
                const labelList& cellFaces = cells[celli];

                forAll(cellFaces, cellFacei)
                {
                    if (patchFace == faces[cellFaces[cellFacei]])
                    {
                        patchFaceOK = true;

                        if
                        (
                            (
                                patchFace.area(points)
                              & faces[cellFaces[cellFacei]].area(points)
                            ) < 0
                        )
                        {
                            oss << "    Face #" << patchFacei << ' '
                                << patchFace << " of patch #" << patchi << ' '
                                << patches[patchi].name() << " points inwards"
                                << endl;
                        }
                    }
                }
            }

            if (!patchFaceOK)
            {
                oss << "    Face #" << patchFacei << ' ' << patchFace
                    << " of patch #" << patchi << ' ' << patches[patchi].name()
                    << " does not match any block faces" << endl;
            }
        }
    }

    if (!oss.str().empty())
    {
        FatalIOErrorInFunction(dict)
            << oss.str().c_str()
            << exit(FatalIOError);
    }
}


// ************************************************************************* //
