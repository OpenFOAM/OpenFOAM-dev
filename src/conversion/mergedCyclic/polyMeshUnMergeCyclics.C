/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2020-2022 OpenFOAM Foundation
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

#include "polyMeshUnMergeCyclics.H"
#include "cyclicPolyPatch.H"
#include "mergedCyclicPolyPatch.H"
#include "polyTopoChange.H"
#include "unitConversion.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

labelList primitivePatchGetZones
(
    const primitivePatch& pp,
    const scalar includedAngle
)
{
    const scalar minCos = Foam::cos(degToRad(180.0 - includedAngle));

    const labelListList& edgeFaces = pp.edgeFaces();
    const labelListList& faceEdges = pp.faceEdges();
    const vectorField& faceNormals = pp.faceNormals();

    // Working arrays
    labelList facePrev(pp.size(), -1);
    labelList faceEdge0(pp.size(), -1);
    labelList faceEdgeFace0(pp.size(), -1);

    // The resulting set of zones
    labelList faceZones(pp.size(), -1);

    // The current zone being populated
    label zonei = 0;

    // The starting face of the current zone
    label facei0 = 0;

    // While there are faces to consider
    while (facei0 < faceZones.size())
    {
        label facei = facei0;
        label faceEdgei = 0;
        label faceEdgeFacei = 0;

        faceZones[facei] = zonei;

        do
        {
            // Find the connected edge and face
            const label edgei = faceEdges[facei][faceEdgei];
            const label facej = edgeFaces[edgei][faceEdgeFacei];

            // Find the corresponding position within the other face
            const label faceEdgej = findIndex(faceEdges[facej], edgei);
            const label faceEdgeFacej = findIndex(edgeFaces[edgei], facej);

            if
            (
                faceZones[facej] == zonei
             || (faceNormals[facei] & faceNormals[facej]) < minCos
             || face::compare(pp[facei], pp[facej]) == -1
             || edge::compare
                (
                    pp[facei].faceEdge(faceEdgei),
                    pp[facej].faceEdge(faceEdgej)
                ) != -1
            )
            {
                // Move to the next face for this edge
                faceEdgeFacei = edgeFaces[edgei].fcIndex(faceEdgeFacei);

                // If wrapped around, move to the next edge on the face
                if (faceEdgeFacei == 0)
                {
                    faceEdgei = faceEdges[facei].fcIndex(faceEdgei);
                }

                // If back at the start then backtrack to the previous face
                if
                (
                    faceEdgei == faceEdge0[facei]
                 && faceEdgeFacei == faceEdgeFace0[facei]
                )
                {
                    facei = facePrev[facei];
                    faceEdgei = faceEdge0[facei];
                    faceEdgeFacei = faceEdgeFace0[facei];
                }
            }
            else
            {
                // Move into the next face
                faceZones[facej] = zonei;
                facePrev[facej] = facei;
                facei = facej;

                // Set the starting position within this face
                faceEdgei = faceEdge0[facej] = faceEdgej;
                faceEdgeFacei = faceEdgeFace0[facej] = faceEdgeFacej;
            }
        }
        while (facei != facei0 || faceEdgei != 0 || faceEdgeFacei != 0);

        // Get the next zone and starting face
        ++ zonei;
        while (facei0 < faceZones.size() && faceZones[facei0] != -1)
        {
            ++ facei0;
        }
    }

    return faceZones;
}


boolList primitivePatchGetHalves
(
    const primitivePatch& pp,
    const scalar includedAngle
)
{
    const labelList faceZones(primitivePatchGetZones(pp, includedAngle));

    // If the number of zones is not two, then the split has failed
    // !!! To improve the robustness of this in the presence of non-planar
    // cyclics we would have to find pairs of zones that all share the same
    // transformation. We would still fail if the number of zones were odd.
    if (findMin(faceZones) != 0 && findMax(faceZones) != 1)
    {
        FatalErrorInFunction
            << "Patch did not divide into halves based on topology and an "
            << "included angle of " << includedAngle << " degrees"
            << exit(FatalError);
    }

    // Convert the zone labels into bools describing which half each face is in
    boolList faceHalves(pp.size());
    forAll(faceHalves, facei)
    {
        faceHalves[facei] = faceZones[facei] == 1;
    }

    return faceHalves;
}

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void Foam::polyMeshUnMergeCyclics(polyMesh& mesh, const scalar includedAngle)
{
    // Add cyclic patches. Clone all the existing patches into a new list, and
    // add two additional empty cyclic patches for each merged-cyclic. Then
    // re-patch the mesh
    DynamicList<polyPatch*> patches;

    forAll(mesh.boundaryMesh(), patchi)
    {
        const polyPatch& pp = mesh.boundaryMesh()[patchi];
        patches.append
        (
            pp.clone
            (
                mesh.boundaryMesh(),
                patches.size(),
                pp.size(),
                pp.start()
            ).ptr()
        );
    }

    bool hasMergedCyclics = false;
    labelList patchHalf0(mesh.boundaryMesh().size(), -1);
    labelList patchHalf1(mesh.boundaryMesh().size(), -1);

    forAll(mesh.boundaryMesh(), patchi)
    {
        const polyPatch& pp = mesh.boundaryMesh()[patchi];
        if (isA<mergedCyclicPolyPatch>(pp))
        {
            hasMergedCyclics = true;

            patchHalf0[patchi] = patches.size();

            patches.append
            (
                new cyclicPolyPatch
                (
                    pp.name() + "-half-0",
                    0,
                    mesh.nFaces(),
                    patches.size(),
                    mesh.boundaryMesh(),
                    cyclicPolyPatch::typeName,
                    pp.name() + "-half-1"
                )
            );

            patchHalf1[patchi] = patches.size();

            patches.append
            (
                new cyclicPolyPatch
                (
                    pp.name() + "-half-1",
                    0,
                    mesh.nFaces(),
                    patches.size(),
                    mesh.boundaryMesh(),
                    cyclicPolyPatch::typeName,
                    pp.name() + "-half-0"
                )
            );
        }
    }

    if (!hasMergedCyclics)
    {
        return;
    }

    mesh.removeBoundary();
    mesh.addPatches(patches);

    // Move faces from the merged cyclics to the pairs of cyclic patches
    polyTopoChange meshMod(mesh, false);
    forAll(mesh.boundaryMesh(), patchi)
    {
        const polyPatch& pp = mesh.boundaryMesh()[patchi];
        if (isA<mergedCyclicPolyPatch>(pp))
        {
            const boolList patchFaceHalves
            (
                primitivePatchGetHalves(pp, includedAngle)
            );

            forAll(pp, patchFacei)
            {
                const label facei = pp.start() + patchFacei;

                meshMod.modifyFace
                (
                    mesh.faces()[facei],
                    facei,
                    mesh.faceOwner()[facei],
                    -1,
                    false,
                    patchFaceHalves[patchFacei]
                  ? patchHalf0[patchi]
                  : patchHalf1[patchi],
                    -1,
                    false
                );
            }
        }
    }
    autoPtr<polyTopoChangeMap> map = meshMod.changeMesh(mesh, false);

    // Remove the (now) empty merged cyclic patches. Copy everything except the
    // merged cyclics into a patch list and use this then re-patch the mesh.
    patches.clear();

    forAll(mesh.boundaryMesh(), patchi)
    {
        const polyPatch& pp = mesh.boundaryMesh()[patchi];
        if (!isA<mergedCyclicPolyPatch>(pp))
        {
            patches.append
            (
                pp.clone
                (
                    mesh.boundaryMesh(),
                    patches.size(),
                    pp.size(),
                    pp.start()
                ).ptr()
            );
        }
    }

    mesh.removeBoundary();
    mesh.addPatches(patches);
}


// ************************************************************************* //
