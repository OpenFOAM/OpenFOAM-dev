/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

// Check the blockMesh topology
void Foam::blockMesh::checkBlockMesh(const polyMesh& bm) const
{
    if (verboseOutput)
    {
        Info<< nl << "Check topology" << endl;
    }

    bool ok = true;

    const pointField& points = bm.points();
    const faceList& faces = bm.faces();
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
                                patchFace.normal(points)
                              & faces[cellFaces[cellFacei]].normal(points)
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
        FatalErrorIn("blockMesh::checkBlockMesh(const polyMesh& bm)")
            << "Block mesh topology incorrect, stopping mesh generation!"
            << exit(FatalError);
    }
}


bool Foam::blockMesh::blockLabelsOK
(
    const label blockLabel,
    const pointField& points,
    const cellShape& blockShape
) const
{
    bool ok = true;

    forAll(blockShape, blockI)
    {
        if (blockShape[blockI] < 0)
        {
            ok = false;

            WarningIn
            (
                "bool Foam::blockMesh::blockLabelsOK(...)"
            )   << "out-of-range point label " << blockShape[blockI]
                << " (min = 0"
                << ") in block " << blockLabel << endl;
        }
        else if (blockShape[blockI] >= points.size())
        {
            ok = false;

            WarningIn
            (
                "bool Foam::blockMesh::blockLabelsOK(...)"
            )   << "out-of-range point label " << blockShape[blockI]
                << " (max = " << points.size() - 1
                << ") in block " << blockLabel << endl;
        }
    }

    return ok;
}


bool Foam::blockMesh::patchLabelsOK
(
    const label patchLabel,
    const pointField& points,
    const faceList& patchFaces
) const
{
    bool ok = true;

    forAll(patchFaces, faceI)
    {
        const labelList& f = patchFaces[faceI];

        forAll(f, fp)
        {
            if (f[fp] < 0)
            {
                ok = false;

                WarningIn
                (
                    "bool Foam::blockMesh::patchLabelsOK(...)"
                )   << "out-of-range point label " << f[fp]
                    << " (min = 0"
                    << ") on patch " << patchLabel
                    << ", face " << faceI << endl;
            }
            else if (f[fp] >= points.size())
            {
                ok = false;

                WarningIn
                (
                    "bool Foam::blockMesh::patchLabelsOK(...)"
                )   << "out-of-range point label " << f[fp]
                    << " (max = " << points.size() - 1
                    << ") on patch " << patchLabel
                    << ", face " << faceI << endl;

            }
        }
    }

    return ok;
}


// ************************************************************************* //
