/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
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

Description
    Test PatchEdgeFaceWave.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "fvMesh.H"
#include "volFields.H"
#include "PatchEdgeFaceWave.H"
#include "patchEdgeFaceInfo.H"
#include "patchPatchDist.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::validArgs.append("patch");

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    // Get name of patch
    const word patchName = args[1];
    const polyPatch& patch = patches[patchName];

    // 1. Walk from a single edge
    {
        // Data on all edges and faces
        List<patchEdgeFaceInfo> allEdgeInfo(patch.nEdges());
        List<patchEdgeFaceInfo> allFaceInfo(patch.size());

        // Initial seed
        DynamicList<label> initialEdges;
        DynamicList<patchEdgeFaceInfo> initialEdgesInfo;


        // Just set an edge on the master
        if (Pstream::master())
        {
            label edgeI = 0;
            Info<< "Starting walk on edge " << edgeI << endl;

            initialEdges.append(edgeI);
            const edge& e = patch.edges()[edgeI];
            initialEdgesInfo.append
            (
                patchEdgeFaceInfo
                (
                    e.centre(patch.localPoints()),
                    0.0
                )
            );
        }


        // Walk
        PatchEdgeFaceWave
        <
            primitivePatch,
            patchEdgeFaceInfo
        > calc
        (
            mesh,
            patch,
            initialEdges,
            initialEdgesInfo,
            allEdgeInfo,
            allFaceInfo,
            returnReduce(patch.nEdges(), sumOp<label>())
        );


        // Extract as patchField
        volScalarField vsf
        (
            IOobject
            (
                "patchDist",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedScalar("patchDist", dimLength, 0.0)
        );
        scalarField pf(vsf.boundaryField()[patch.index()].size());
        forAll(pf, facei)
        {
            pf[facei] = Foam::sqrt(allFaceInfo[facei].distSqr());
        }
        vsf.boundaryFieldRef()[patch.index()] = pf;

        Info<< "Writing patchDist volScalarField to " << runTime.value()
            << endl;

        vsf.write();
    }


    // 2. Use a wrapper to walk from all boundary edges on selected patches
    {
        labelHashSet otherPatchIDs(identity(mesh.boundaryMesh().size()));
        otherPatchIDs.erase(patch.index());

        Info<< "Walking on patch " << patch.index()
            << " from edges shared with patches " << otherPatchIDs
            << endl;

        patchPatchDist pwd(patch, otherPatchIDs);

        // Extract as patchField
        volScalarField vsf
        (
            IOobject
            (
                "otherPatchDist",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedScalar("otherPatchDist", dimLength, 0.0)
        );
        vsf.boundaryFieldRef()[patch.index()] = pwd;

        Info<< "Writing otherPatchDist volScalarField to " << runTime.value()
            << endl;

        vsf.write();
    }

    Info<< "\nEnd\n" << endl;
    return 0;
}


// ************************************************************************* //
