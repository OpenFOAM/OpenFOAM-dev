/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013 OpenFOAM Foundation
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

Description
    Detect point pinches

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "PatchTools.H"
#include "Time.H"
#include "polyMesh.H"
#include "patchEdgeFaceRegions.H"
#include "PatchEdgeFaceWave.H"
#include "globalIndex.H"
#include "syncTools.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// Main program:

int main(int argc, char *argv[])
{
   argList::validArgs.append("patch");

#   include "setRootCase.H"
#   include "createTime.H"

    const word patchName = args[1];

#   include "createPolyMesh.H"

    Info<< "Mesh read in = "
        << runTime.cpuTimeIncrement()
        << " s\n" << endl << endl;


    const polyBoundaryMesh& pbm = mesh.boundaryMesh();
    label patchI = pbm.findPatchID(patchName);
    const polyPatch& patch = pbm[patchI];

    Info<< "Patch:" << patch.name() << endl;



    // Data on all edges and faces
    List<patchEdgeFaceRegions> allEdgeInfo(patch.nEdges());
    List<patchEdgeFaceRegions> allFaceInfo(patch.size());

    // Determine parallel global indexing
    const globalIndex globalNumbering(patch.size());

    DynamicList<label> changedEdges(4*patch.size());
    DynamicList<patchEdgeFaceRegions> changedInfo(changedEdges.size());

    forAll(patch, faceI)
    {
        const labelList& fEdges = patch.faceEdges()[faceI];

        label globalFaceI = globalNumbering.toGlobal(faceI);

        forAll(fEdges, i)
        {
            changedEdges.append(fEdges[i]);
            changedInfo.append
            (
                patchEdgeFaceRegions(labelPair(globalFaceI, globalFaceI))
            );
        }
    }

    // Walk
    PatchEdgeFaceWave
    <
        primitivePatch,
        patchEdgeFaceRegions
    > calc
    (
        mesh,
        patch,
        changedEdges,
        changedInfo,
        allEdgeInfo,
        allFaceInfo,
        returnReduce(patch.nEdges(), sumOp<label>())
    );


    Info<< "Time now = " << runTime.timeName() << endl;


    // Detect points with multiple regions
    labelList duplicateRegion(patch.nPoints(), -1);
    {
        labelList currentRegion(patch.nPoints(), -1);

        forAll(patch.localFaces(), faceI)
        {
            const face& f = patch.localFaces()[faceI];

            forAll(f, fp)
            {
                label faceRegion = allFaceInfo[faceI].regions()[fp];

                label pointI = f[fp];

                if (currentRegion[pointI] == -1)
                {
                    currentRegion[pointI] = faceRegion;
                }
                else if (currentRegion[pointI] != faceRegion)
                {
                    if (duplicateRegion[pointI] == -1)
                    {
                        Pout<< "Multi region point:"
                            << patch.localPoints()[pointI]
                            << " with region:" << currentRegion[pointI]
                            << " and region:" << faceRegion
                            << endl;
                        duplicateRegion[pointI] = currentRegion[pointI];
                    }
                }
            }
        }
    }


    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
