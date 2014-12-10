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

#include "polyMeshTools.H"
#include "syncTools.H"
#include "pyramidPointFaceRef.H"
#include "primitiveMeshTools.H"
#include "polyMeshTools.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::scalarField> Foam::polyMeshTools::faceOrthogonality
(
    const polyMesh& mesh,
    const vectorField& areas,
    const vectorField& cc
)
{
    const labelList& own = mesh.faceOwner();
    const labelList& nei = mesh.faceNeighbour();
    const polyBoundaryMesh& pbm = mesh.boundaryMesh();

    tmp<scalarField> tortho(new scalarField(mesh.nFaces(), 1.0));
    scalarField& ortho = tortho();

    // Internal faces
    forAll(nei, faceI)
    {
        ortho[faceI] = primitiveMeshTools::faceOrthogonality
        (
            cc[own[faceI]],
            cc[nei[faceI]],
            areas[faceI]
        );
    }


    // Coupled faces

    pointField neighbourCc;
    syncTools::swapBoundaryCellPositions(mesh, cc, neighbourCc);

    forAll(pbm, patchI)
    {
        const polyPatch& pp = pbm[patchI];
        if (pp.coupled())
        {
            forAll(pp, i)
            {
                label faceI = pp.start() + i;
                label bFaceI = faceI - mesh.nInternalFaces();

                ortho[faceI] = primitiveMeshTools::faceOrthogonality
                (
                    cc[own[faceI]],
                    neighbourCc[bFaceI],
                    areas[faceI]
                );
            }
        }
    }

    return tortho;
}


Foam::tmp<Foam::scalarField> Foam::polyMeshTools::faceSkewness
(
    const polyMesh& mesh,
    const pointField& p,
    const vectorField& fCtrs,
    const vectorField& fAreas,
    const vectorField& cellCtrs
)
{
    const labelList& own = mesh.faceOwner();
    const labelList& nei = mesh.faceNeighbour();
    const polyBoundaryMesh& pbm = mesh.boundaryMesh();

    tmp<scalarField> tskew(new scalarField(mesh.nFaces()));
    scalarField& skew = tskew();

    forAll(nei, faceI)
    {
        skew[faceI] = primitiveMeshTools::faceSkewness
        (
            mesh,
            p,
            fCtrs,
            fAreas,

            faceI,
            cellCtrs[own[faceI]],
            cellCtrs[nei[faceI]]
        );
    }


    // Boundary faces: consider them to have only skewness error.
    // (i.e. treat as if mirror cell on other side)

    pointField neighbourCc;
    syncTools::swapBoundaryCellPositions(mesh, cellCtrs, neighbourCc);

    forAll(pbm, patchI)
    {
        const polyPatch& pp = pbm[patchI];
        if (pp.coupled())
        {
            forAll(pp, i)
            {
                label faceI = pp.start() + i;
                label bFaceI = faceI - mesh.nInternalFaces();

                skew[faceI] = primitiveMeshTools::faceSkewness
                (
                    mesh,
                    p,
                    fCtrs,
                    fAreas,

                    faceI,
                    cellCtrs[own[faceI]],
                    neighbourCc[bFaceI]
                );
            }
        }
        else
        {
            forAll(pp, i)
            {
                label faceI = pp.start() + i;

                skew[faceI] = primitiveMeshTools::boundaryFaceSkewness
                (
                    mesh,
                    p,
                    fCtrs,
                    fAreas,

                    faceI,
                    cellCtrs[own[faceI]]
                );
            }
        }
    }

    return tskew;
}


Foam::tmp<Foam::scalarField> Foam::polyMeshTools::faceWeights
(
    const polyMesh& mesh,
    const vectorField& fCtrs,
    const vectorField& fAreas,
    const vectorField& cellCtrs
)
{
    const labelList& own = mesh.faceOwner();
    const labelList& nei = mesh.faceNeighbour();
    const polyBoundaryMesh& pbm = mesh.boundaryMesh();

    tmp<scalarField> tweight(new scalarField(mesh.nFaces(), 1.0));
    scalarField& weight = tweight();

    // Internal faces
    forAll(nei, faceI)
    {
        const point& fc = fCtrs[faceI];
        const vector& fa = fAreas[faceI];

        scalar dOwn = mag(fa & (fc-cellCtrs[own[faceI]]));
        scalar dNei = mag(fa & (cellCtrs[nei[faceI]]-fc));

        weight[faceI] = min(dNei,dOwn)/(dNei+dOwn+VSMALL);
    }


    // Coupled faces

    pointField neiCc;
    syncTools::swapBoundaryCellPositions(mesh, cellCtrs, neiCc);

    forAll(pbm, patchI)
    {
        const polyPatch& pp = pbm[patchI];
        if (pp.coupled())
        {
            forAll(pp, i)
            {
                label faceI = pp.start() + i;
                label bFaceI = faceI - mesh.nInternalFaces();

                const point& fc = fCtrs[faceI];
                const vector& fa = fAreas[faceI];

                scalar dOwn = mag(fa & (fc-cellCtrs[own[faceI]]));
                scalar dNei = mag(fa & (neiCc[bFaceI]-fc));

                weight[faceI] = min(dNei,dOwn)/(dNei+dOwn+VSMALL);
            }
        }
    }

    return tweight;
}


Foam::tmp<Foam::scalarField> Foam::polyMeshTools::volRatio
(
    const polyMesh& mesh,
    const scalarField& vol
)
{
    const labelList& own = mesh.faceOwner();
    const labelList& nei = mesh.faceNeighbour();
    const polyBoundaryMesh& pbm = mesh.boundaryMesh();

    tmp<scalarField> tratio(new scalarField(mesh.nFaces(), 1.0));
    scalarField& ratio = tratio();

    // Internal faces
    forAll(nei, faceI)
    {
        scalar volOwn = vol[own[faceI]];
        scalar volNei = vol[nei[faceI]];

        ratio[faceI] = min(volOwn,volNei)/(max(volOwn, volNei)+VSMALL);
    }


    // Coupled faces

    scalarField neiVol;
    syncTools::swapBoundaryCellList(mesh, vol, neiVol);

    forAll(pbm, patchI)
    {
        const polyPatch& pp = pbm[patchI];
        if (pp.coupled())
        {
            forAll(pp, i)
            {
                label faceI = pp.start() + i;
                label bFaceI = faceI - mesh.nInternalFaces();

                scalar volOwn = vol[own[faceI]];
                scalar volNei = neiVol[bFaceI];

                ratio[faceI] = min(volOwn,volNei)/(max(volOwn, volNei)+VSMALL);
            }
        }
    }

    return tratio;
}


// ************************************************************************* //
