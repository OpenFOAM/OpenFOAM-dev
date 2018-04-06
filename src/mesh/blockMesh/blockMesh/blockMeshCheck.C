/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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
    Info<< nl << "Check topology" << endl;

    bool ok = true;

    // Check for duplicate curved edge definitions
    forAll(edges_, cei)
    {
        for (label cej=cei+1; cej<edges_.size(); cej++)
        {
            if (edges_[cei].compare(edges_[cej]) != 0)
            {
                Info<< "    Curved edge ";
                edges_[cej].write(Info, dict);
                Info<< "    is a duplicate of curved edge " << edges_[cei]
                    << endl;
                ok = false;
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
            Info<< "    Curved edge ";
            edges_[cei].write(Info, dict);
            Info<< "    does not correspond to a block edge."
                << endl;
            ok = false;
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
                Info<< "    Curved face ";
                faces_[cfj].write(Info, dict);
                Info<< "    is a duplicate of curved face ";
                faces_[cfi].write(Info, dict);
                Info<< endl;
                ok = false;
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
            Info<< "    Curved face ";
            faces_[cfi].write(Info, dict);
            Info<< "    does not correspond to a block face." << endl;
            ok = false;
        }
    }

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
        Info<< nl << tab << "Basic statistics" << nl
            << tab << tab << "Number of internal faces : "
            << bm.nInternalFaces() << nl
            << tab << tab << "Number of boundary faces : "
            << nBoundaryFaces << nl
            << tab << tab << "Number of defined boundary faces : "
            << nDefinedBoundaryFaces << nl
            << tab << tab << "Number of undefined boundary faces : "
            << nBoundaryFaces - nDefinedBoundaryFaces << nl;

        if ((nBoundaryFaces - nDefinedBoundaryFaces) > 0)
        {
            Info<< tab << tab << tab
                << "(Warning : only leave undefined the front and back planes "
                << "of 2D planar geometries!)" << endl;
        }

        Info<< tab << "Checking patch -> block consistency" << endl;
    }


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
                            ) < 0.0
                        )
                        {
                            Info<< tab << tab
                                << "Face " << patchFacei
                                << " of patch " << patchi
                                << " (" << patches[patchi].name() << ")"
                                << " points inwards"
                                << endl;

                            ok = false;
                        }
                    }
                }
            }

            if (!patchFaceOK)
            {
                Info<< tab << tab
                    << "Face " << patchFacei
                    << " of patch " << patchi
                    << " (" << patches[patchi].name() << ")"
                    << " does not match any block faces" << endl;

                ok = false;
            }
        }
    }

    if (verboseOutput)
    {
        Info<< endl;
    }

    if (!ok)
    {
        FatalErrorInFunction
            << "Block mesh topology incorrect, stopping mesh generation!"
            << exit(FatalError);
    }
}


// ************************************************************************* //
