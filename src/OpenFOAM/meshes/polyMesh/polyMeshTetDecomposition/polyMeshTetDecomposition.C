/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2023 OpenFOAM Foundation
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

#include "polyMeshTetDecomposition.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

// Note: the use of this tolerance is ad-hoc, there may be extreme
// cases where the resulting tetrahedra still have particle tracking
// problems, or tets with lower quality may track OK.
const Foam::scalar Foam::polyMeshTetDecomposition::minTetQuality = sqr(small);


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

Foam::scalar Foam::polyMeshTetDecomposition::minQuality
(
    const polyMesh& mesh,
    const point& cC,
    const label fI,
    const bool isOwner,
    const label faceBasePtI
)
{
    // Does fan decomposition of face (starting at faceBasePti) and determines
    // min quality over all resulting tets.

    const pointField& pPts = mesh.points();
    const face& f = mesh.faces()[fI];
    const point& tetBasePt = pPts[f[faceBasePtI]];

    scalar thisBaseMinTetQuality = vGreat;

    for (label tetPtI = 1; tetPtI < f.size() - 1; tetPtI++)
    {
        label facePtI = (tetPtI + faceBasePtI) % f.size();
        label otherFacePtI = f.fcIndex(facePtI);

        label ptAI = -1;
        label ptBI = -1;

        if (isOwner)
        {
            ptAI = f[facePtI];
            ptBI = f[otherFacePtI];
        }
        else
        {
            ptAI = f[otherFacePtI];
            ptBI = f[facePtI];
        }

        const point& pA = pPts[ptAI];
        const point& pB = pPts[ptBI];

        tetPointRef tet(cC, tetBasePt, pA, pB);

        scalar tetQuality = tet.quality();

        if (tetQuality < thisBaseMinTetQuality)
        {
            thisBaseMinTetQuality = tetQuality;
        }
    }
    return thisBaseMinTetQuality;
}


Foam::label Foam::polyMeshTetDecomposition::findSharedBasePoint
(
    const polyMesh& mesh,
    label fI,
    const point& nCc,
    scalar tol,
    bool report
)
{
    const faceList& pFaces = mesh.faces();
    const vectorField& pC = mesh.cellCentres();
    const labelList& pOwner = mesh.faceOwner();

    const face& f = pFaces[fI];

    label oCI = pOwner[fI];

    const point& oCc = pC[oCI];

    forAll(f, faceBasePtI)
    {
        scalar minQ = minQuality(mesh, oCc, fI, true, faceBasePtI);
        minQ = min(minQ, minQuality(mesh, nCc, fI, false, faceBasePtI));

        if (minQ > tol)
        {
            return faceBasePtI;
        }
    }

    // If a base point hasn't triggered a return by now, then there is
    // none that can produce a good decomposition
    return -1;
}


Foam::label Foam::polyMeshTetDecomposition::findSharedBasePoint
(
    const polyMesh& mesh,
    label fI,
    scalar tol,
    bool report
)
{
    return findSharedBasePoint
    (
        mesh,
        fI,
        mesh.cellCentres()[mesh.faceNeighbour()[fI]],
        tol,
        report
    );
}


Foam::label Foam::polyMeshTetDecomposition::findBasePoint
(
    const polyMesh& mesh,
    label fI,
    scalar tol,
    bool report
)
{
    const faceList& pFaces = mesh.faces();
    const vectorField& pC = mesh.cellCentres();
    const labelList& pOwner = mesh.faceOwner();

    const face& f = pFaces[fI];

    label cI = pOwner[fI];

    bool own = (pOwner[fI] == cI);

    const point& cC = pC[cI];

    forAll(f, faceBasePtI)
    {
        scalar quality = minQuality(mesh, cC, fI, own, faceBasePtI);

        if (quality > tol)
        {
            return faceBasePtI;
        }
    }

    // If a base point hasn't triggered a return by now, then there is
    // none that can produce a good decomposition
    return -1;
}


Foam::labelList Foam::polyMeshTetDecomposition::findFaceBasePts
(
    const polyMesh& mesh,
    scalar tol,
    bool report
)
{
    const labelList& pOwner = mesh.faceOwner();
    const vectorField& pC = mesh.cellCentres();

    // Find a suitable base point for each face, considering both
    // cells for interface faces or those on coupled patches

    labelList tetBasePtIs(mesh.nFaces(), -1);

    label nInternalFaces = mesh.nInternalFaces();

    for (label fI = 0; fI < nInternalFaces; fI++)
    {
        tetBasePtIs[fI] = findSharedBasePoint(mesh, fI, tol, report);
    }

    pointField neighbourCellCentres(mesh.nFaces() - nInternalFaces);

    for(label facei = nInternalFaces; facei < mesh.nFaces(); facei++)
    {
        neighbourCellCentres[facei - nInternalFaces] =
            pC[pOwner[facei]];
    }

    syncTools::swapBoundaryFacePositions(mesh, neighbourCellCentres);

    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    SubList<label> boundaryFaceTetBasePtIs
    (
        tetBasePtIs,
        mesh.nFaces() - nInternalFaces,
        nInternalFaces
    );

    for
    (
        label fI = nInternalFaces, bFI = 0;
        fI < mesh.nFaces();
        fI++, bFI++
    )
    {
        const label patchi =
            mesh.boundaryMesh().patchIndices()[bFI];

        if (patches[patchi].coupled())
        {
            const coupledPolyPatch& pp =
                refCast<const coupledPolyPatch>(patches[patchi]);

            if (pp.owner())
            {
                boundaryFaceTetBasePtIs[bFI] = findSharedBasePoint
                (
                    mesh,
                    fI,
                    neighbourCellCentres[bFI],
                    tol,
                    report
                );
            }
            else
            {
                // Assign -2, to distinguish from a failed base point
                // find, which returns -1.
                boundaryFaceTetBasePtIs[bFI] = -2;
            }
        }
        else
        {
            boundaryFaceTetBasePtIs[bFI] = findBasePoint
            (
                mesh,
                fI,
                tol,
                report
            );
        }
    }

    // maxEqOp will replace the -2 values on the neighbour patches
    // with the result from the owner base point find.

    syncTools::syncBoundaryFaceList
    (
        mesh,
        boundaryFaceTetBasePtIs,
        maxEqOp<label>()
    );

    for
    (
        label fI = nInternalFaces, bFI = 0;
        fI < mesh.nFaces();
        fI++, bFI++
    )
    {
        label& bFTetBasePtI = boundaryFaceTetBasePtIs[bFI];

        if (bFTetBasePtI == -2)
        {
            FatalErrorInFunction
                << "Coupled face base point exchange failure for face "
                << fI
                << abort(FatalError);
        }

        if (bFTetBasePtI < 1)
        {
            // If the base point is -1, it should be left as such to
            // indicate a problem, if it is 0, then no action is required.

            continue;
        }

        const label patchi = mesh.boundaryMesh().patchIndices()[bFI];

        if (patches[patchi].coupled())
        {
            const coupledPolyPatch& pp =
                refCast<const coupledPolyPatch>(patches[patchi]);

            // Calculated base points on coupled faces are those of
            // the owner patch face. They need to be reindexed to for
            // the non-owner face, which has the opposite order.

            // So, for fPtI_o != 0, fPtI_n = f.size() - fPtI_o

            // i.e.:

            // owner coupledPolyPatch face
            // face    (a b c d e f)
            // fPtI     0 1 2 3 4 5
            //            +
            //              #

            // neighbour coupledPolyPatch face
            // face    (a f e d c b)
            // fPtI     0 1 2 3 4 5
            //                    +
            //                  #
            // +: 6 - 1 = 5
            // #: 6 - 2 = 4

            if (!pp.owner())
            {
                bFTetBasePtI = mesh.faces()[fI].size() - bFTetBasePtI;
            }
        }
    }

    return tetBasePtIs;
}


bool Foam::polyMeshTetDecomposition::checkFaceTets
(
    const polyMesh& mesh,
    scalar tol,
    const bool report,
    labelHashSet* setPtr
)
{
    const labelList& own = mesh.faceOwner();
    const labelList& nei = mesh.faceNeighbour();
    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    const vectorField& cc = mesh.cellCentres();
    const vectorField& fc = mesh.faceCentres();

    // Calculate coupled cell centre
    pointField neiCc(mesh.nFaces() - mesh.nInternalFaces());

    for (label facei = mesh.nInternalFaces(); facei < mesh.nFaces(); facei++)
    {
        neiCc[facei - mesh.nInternalFaces()] = cc[own[facei]];
    }

    syncTools::swapBoundaryFacePositions(mesh, neiCc);

    const faceList& fcs = mesh.faces();

    const pointField& p = mesh.points();

    label nErrorTets = 0;

    forAll(fcs, facei)
    {
        const face& f = fcs[facei];

        forAll(f, fPtI)
        {
            scalar tetQual = tetPointRef
            (
                p[f[fPtI]],
                p[f.nextLabel(fPtI)],
                fc[facei],
                cc[own[facei]]
            ).quality();

            if (tetQual > -tol)
            {
                if (setPtr)
                {
                    setPtr->insert(facei);
                }

                nErrorTets++;
                break;              // no need to check other tets
            }
        }

        if (mesh.isInternalFace(facei))
        {
            // Create the neighbour tet - it will have positive volume
            const face& f = fcs[facei];

            forAll(f, fPtI)
            {
                scalar tetQual = tetPointRef
                (
                    p[f[fPtI]],
                    p[f.nextLabel(fPtI)],
                    fc[facei],
                    cc[nei[facei]]
                ).quality();

                if (tetQual < tol)
                {
                    if (setPtr)
                    {
                        setPtr->insert(facei);
                    }

                    nErrorTets++;
                    break;
                }
            }

            if (findSharedBasePoint(mesh, facei, tol, report) == -1)
            {
                if (setPtr)
                {
                    setPtr->insert(facei);
                }

                nErrorTets++;
            }
        }
        else
        {
            const label patchi =
                patches.patchIndices()[facei - mesh.nInternalFaces()];

            if (patches[patchi].coupled())
            {
                if
                (
                    findSharedBasePoint
                    (
                        mesh,
                        facei,
                        neiCc[facei - mesh.nInternalFaces()],
                        tol,
                        report
                    ) == -1
                )
                {
                    if (setPtr)
                    {
                        setPtr->insert(facei);
                    }

                    nErrorTets++;
                }
            }
            else
            {
                if (findBasePoint(mesh, facei, tol, report) == -1)
                {
                    if (setPtr)
                    {
                        setPtr->insert(facei);
                    }

                    nErrorTets++;
                }
            }
        }
    }

    reduce(nErrorTets, sumOp<label>());

    if (nErrorTets > 0)
    {
        if (report)
        {
            Info<< " ***Error in face tets: "
                << nErrorTets << " faces with low quality or negative volume "
                << "decomposition tets." << endl;
        }

        return true;
    }
    else
    {
        if (report)
        {
            Info<< "    Face tets OK." << endl;
        }

        return false;
    }
}


Foam::List<Foam::tetIndices> Foam::polyMeshTetDecomposition::faceTetIndices
(
    const polyMesh& mesh,
    label fI,
    label cI
)
{
    const faceList& pFaces = mesh.faces();

    const face& f = pFaces[fI];

    label nTets = f.size() - 2;

    List<tetIndices> faceTets(nTets);

    for (label tetPtI = 1; tetPtI < f.size() - 1; tetPtI ++)
    {
        faceTets[tetPtI - 1] = tetIndices(cI, fI, tetPtI);
    }

    return faceTets;
}


Foam::List<Foam::tetIndices> Foam::polyMeshTetDecomposition::cellTetIndices
(
    const polyMesh& mesh,
    label cI
)
{
    const faceList& pFaces = mesh.faces();
    const cellList& pCells = mesh.cells();

    const cell& thisCell = pCells[cI];

    label nTets = 0;

    forAll(thisCell, cFI)
    {
        nTets += pFaces[thisCell[cFI]].size() - 2;
    }

    DynamicList<tetIndices> cellTets(nTets);

    forAll(thisCell, cFI)
    {
        label fI = thisCell[cFI];

        cellTets.append(faceTetIndices(mesh, fI, cI));
    }

    return cellTets;
}


Foam::tetIndices Foam::polyMeshTetDecomposition::findTet
(
    const polyMesh& mesh,
    label cI,
    const point& pt
)
{
    const faceList& pFaces = mesh.faces();
    const cellList& pCells = mesh.cells();

    const cell& thisCell = pCells[cI];

    tetIndices tetContainingPt;


    forAll(thisCell, cFI)
    {
        label fI = thisCell[cFI];
        const face& f = pFaces[fI];

        for (label tetPtI = 1; tetPtI < f.size() - 1; tetPtI++)
        {
            // Get tetIndices of face triangle
            tetIndices faceTetIs(cI, fI, tetPtI);

            // Check if inside
            if (faceTetIs.tet(mesh).inside(pt))
            {
                tetContainingPt = faceTetIs;
                break;
            }
        }

        if (tetContainingPt.cell() != -1)
        {
            break;
        }
    }

    return tetContainingPt;
}


// ************************************************************************* //
