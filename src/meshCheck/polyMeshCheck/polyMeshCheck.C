/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2012-2023 OpenFOAM Foundation
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

#include "primitiveMeshCheck.H"
#include "polyMeshCheck.H"
#include "unitConversion.H"
#include "syncTools.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::scalarField> Foam::meshCheck::faceOrthogonality
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
    scalarField& ortho = tortho.ref();

    // Internal faces
    forAll(nei, facei)
    {
        ortho[facei] = meshCheck::faceOrthogonality
        (
            cc[own[facei]],
            cc[nei[facei]],
            areas[facei]
        );
    }


    // Coupled faces

    pointField neighbourCc;
    syncTools::swapBoundaryCellPositions(mesh, cc, neighbourCc);

    forAll(pbm, patchi)
    {
        const polyPatch& pp = pbm[patchi];
        if (pp.coupled())
        {
            forAll(pp, i)
            {
                label facei = pp.start() + i;
                label bFacei = facei - mesh.nInternalFaces();

                ortho[facei] = meshCheck::faceOrthogonality
                (
                    cc[own[facei]],
                    neighbourCc[bFacei],
                    areas[facei]
                );
            }
        }
    }

    return tortho;
}


Foam::tmp<Foam::scalarField> Foam::meshCheck::faceSkewness
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
    scalarField& skew = tskew.ref();

    forAll(nei, facei)
    {
        skew[facei] = meshCheck::faceSkewness
        (
            mesh,
            p,
            fCtrs,
            fAreas,

            facei,
            cellCtrs[own[facei]],
            cellCtrs[nei[facei]]
        );
    }


    // Boundary faces: consider them to have only skewness error.
    // (i.e. treat as if mirror cell on other side)

    pointField neighbourCc;
    syncTools::swapBoundaryCellPositions(mesh, cellCtrs, neighbourCc);

    forAll(pbm, patchi)
    {
        const polyPatch& pp = pbm[patchi];
        if (pp.coupled())
        {
            forAll(pp, i)
            {
                const label facei = pp.start() + i;
                const label bFacei = facei - mesh.nInternalFaces();

                skew[facei] = meshCheck::faceSkewness
                (
                    mesh,
                    p,
                    fCtrs,
                    fAreas,

                    facei,
                    cellCtrs[own[facei]],
                    neighbourCc[bFacei]
                );
            }
        }
        else
        {
            forAll(pp, i)
            {
                const label facei = pp.start() + i;

                skew[facei] = meshCheck::boundaryFaceSkewness
                (
                    mesh,
                    p,
                    fCtrs,
                    fAreas,

                    facei,
                    cellCtrs[own[facei]]
                );
            }
        }
    }

    return tskew;
}


Foam::tmp<Foam::scalarField> Foam::meshCheck::faceWeights
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
    scalarField& weight = tweight.ref();

    // Internal faces
    forAll(nei, facei)
    {
        const point& fc = fCtrs[facei];
        const vector& fa = fAreas[facei];

        const scalar dOwn = mag(fa & (fc-cellCtrs[own[facei]]));
        const scalar dNei = mag(fa & (cellCtrs[nei[facei]]-fc));

        weight[facei] = min(dNei,dOwn)/(dNei+dOwn+vSmall);
    }


    // Coupled faces

    pointField neiCc;
    syncTools::swapBoundaryCellPositions(mesh, cellCtrs, neiCc);

    forAll(pbm, patchi)
    {
        const polyPatch& pp = pbm[patchi];
        if (pp.coupled())
        {
            forAll(pp, i)
            {
                const label facei = pp.start() + i;
                const label bFacei = facei - mesh.nInternalFaces();

                const point& fc = fCtrs[facei];
                const vector& fa = fAreas[facei];

                const scalar dOwn = mag(fa & (fc-cellCtrs[own[facei]]));
                const scalar dNei = mag(fa & (neiCc[bFacei]-fc));

                weight[facei] = min(dNei,dOwn)/(dNei+dOwn+vSmall);
            }
        }
    }

    return tweight;
}


Foam::tmp<Foam::scalarField> Foam::meshCheck::volRatio
(
    const polyMesh& mesh,
    const scalarField& vol
)
{
    const labelList& own = mesh.faceOwner();
    const labelList& nei = mesh.faceNeighbour();
    const polyBoundaryMesh& pbm = mesh.boundaryMesh();

    tmp<scalarField> tratio(new scalarField(mesh.nFaces(), 1.0));
    scalarField& ratio = tratio.ref();

    // Internal faces
    forAll(nei, facei)
    {
        const scalar volOwn = vol[own[facei]];
        const scalar volNei = vol[nei[facei]];

        ratio[facei] = min(volOwn,volNei)/(max(volOwn, volNei)+vSmall);
    }


    // Coupled faces

    scalarField neiVol;
    syncTools::swapBoundaryCellList(mesh, vol, neiVol);

    forAll(pbm, patchi)
    {
        const polyPatch& pp = pbm[patchi];
        if (pp.coupled())
        {
            forAll(pp, i)
            {
                const label facei = pp.start() + i;
                const label bFacei = facei - mesh.nInternalFaces();

                const scalar volOwn = vol[own[facei]];
                const scalar volNei = neiVol[bFacei];

                ratio[facei] = min(volOwn,volNei)/(max(volOwn, volNei)+vSmall);
            }
        }
    }

    return tratio;
}


bool Foam::meshCheck::checkFaceOrthogonality
(
    const polyMesh& mesh,
    const scalar nonOrthThreshold,
    const bool report,
    labelHashSet* setPtr
)
{
    if (mesh.debug)
    {
        InfoInFunction << "Checking mesh non-orthogonality" << endl;
    }

    const vectorField& fAreas = mesh.faceAreas();
    const vectorField& cellCtrs = mesh.cellCentres();

    // Calculate orthogonality for all internal and coupled boundary faces
    // (1 for uncoupled boundary faces)
    tmp<scalarField> tortho = meshCheck::faceOrthogonality
    (
        mesh,
        fAreas,
        cellCtrs
    );
    const scalarField& ortho = tortho.ref();

    // Severe nonorthogonality threshold
    const scalar severeNonorthogonalityThreshold =
        ::cos(degToRad(nonOrthThreshold));


    scalar minDDotS = great;
    scalar sumDDotS = 0.0;
    label nSummed = 0;
    label severeNonOrth = 0;
    label errorNonOrth = 0;


    // Statistics only for internal and masters of coupled faces
    PackedBoolList isMasterFace(syncTools::getInternalOrMasterFaces(mesh));

    forAll(ortho, facei)
    {
        if (ortho[facei] < severeNonorthogonalityThreshold)
        {
            if (ortho[facei] > small)
            {
                if (setPtr)
                {
                    setPtr->insert(facei);
                }

                severeNonOrth++;
            }
            else
            {
                // Error : non-ortho too large
                if (setPtr)
                {
                    setPtr->insert(facei);
                }

                errorNonOrth++;
            }
        }

        if (isMasterFace[facei])
        {
            minDDotS = min(minDDotS, ortho[facei]);
            sumDDotS += ortho[facei];
            nSummed++;
        }
    }

    reduce(minDDotS, minOp<scalar>());
    reduce(sumDDotS, sumOp<scalar>());
    reduce(nSummed, sumOp<label>());
    reduce(severeNonOrth, sumOp<label>());
    reduce(errorNonOrth, sumOp<label>());

    if (report)
    {
        if (nSummed > 0)
        {
            if (report)
            {
                Info<< "    Mesh non-orthogonality Max: "
                    << radToDeg(::acos(min(1.0, max(-1.0, minDDotS))))
                    << " average: "
                    << radToDeg(::acos(min(1.0, max(-1.0, sumDDotS/nSummed))))
                    << endl;
            }
        }

        if (severeNonOrth > 0)
        {
            Info<< "   *Number of severely non-orthogonal (> "
                << nonOrthThreshold << " degrees) faces: "
                << severeNonOrth << "." << endl;
        }
    }

    if (errorNonOrth > 0)
    {
        if (report)
        {
            Info<< " ***Number of non-orthogonality errors: "
                << errorNonOrth << "." << endl;
        }

        return true;
    }
    else
    {
        if (report)
        {
            Info<< "    Non-orthogonality check OK." << endl;
        }

        return false;
    }
}


bool Foam::meshCheck::checkFaceSkewness
(
    const polyMesh& mesh,
    const scalar skewThreshold,
    const bool report,
    labelHashSet* setPtr
)
{
    if (mesh.debug)
    {
        InfoInFunction << "Checking face skewness" << endl;
    }

    const pointField& points = mesh.points();
    const vectorField& fCtrs = mesh.faceCentres();
    const vectorField& fAreas = mesh.faceAreas();
    const vectorField& cellCtrs = mesh.cellCentres();

    // Warn if the skew correction vector is more than skewWarning times
    // larger than the face area vector

    tmp<scalarField> tskew = meshCheck::faceSkewness
    (
        mesh,
        points,
        fCtrs,
        fAreas,
        cellCtrs
    );
    const scalarField& skew = tskew.ref();

    scalar maxSkew = max(skew);
    label nWarnSkew = 0;

    // Statistics only for all faces except slave coupled faces
    const PackedBoolList isMasterFace(syncTools::getMasterFaces(mesh));

    forAll(skew, facei)
    {
        // Check if the skewness vector is greater than the PN vector.
        // This does not cause trouble but is a good indication of a poor mesh.
        if (skew[facei] > skewThreshold)
        {
            if (setPtr)
            {
                setPtr->insert(facei);
            }

            if (isMasterFace[facei])
            {
                nWarnSkew++;
            }
        }
    }

    reduce(maxSkew, maxOp<scalar>());
    reduce(nWarnSkew, sumOp<label>());

    if (nWarnSkew > 0)
    {
        if (report)
        {
            Info<< " ***Max skewness = " << maxSkew
                << ", " << nWarnSkew << " highly skew faces detected"
                   " which may impair the quality of the results"
                << endl;
        }

        return true;
    }
    else
    {
        if (report)
        {
            Info<< "    Max skewness = " << maxSkew << " OK." << endl;
        }

        return false;
    }
}


bool Foam::meshCheck::checkEdgeAlignment
(
    const polyMesh& mesh,
    const bool report,
    const Vector<label>& directions,
    labelHashSet* setPtr
)
{
    // Check 1D/2Dness of edges. Gets passed the non-empty directions and
    // checks all edges in the mesh whether they:
    // - have no component in a non-empty direction or
    // - are only in a single non-empty direction.
    // Empty direction info is passed in as a vector of labels (synchronised)
    // which are 1 if the direction is non-empty, 0 if it is.

    if (mesh.debug)
    {
        InfoInFunction << "Checking edge alignment" << endl;
    }

    const pointField& p = mesh.points();

    label nDirs = 0;
    for (direction cmpt=0; cmpt<vector::nComponents; cmpt++)
    {
        if (directions[cmpt] == 1)
        {
            nDirs++;
        }
        else if (directions[cmpt] != 0)
        {
            FatalErrorInFunction
                << "directions should contain 0 or 1 but is now " << directions
                << exit(FatalError);
        }
    }

    if (nDirs == vector::nComponents)
    {
        return false;
    }


    const faceList& fcs = mesh.faces();

    EdgeMap<label> edgesInError;

    forAll(fcs, facei)
    {
        const face& f = fcs[facei];

        forAll(f, fp)
        {
            const label p0 = f[fp];
            const label p1 = f.nextLabel(fp);

            if (p0 < p1)
            {
                vector d(p[p1]-p[p0]);
                const scalar magD = mag(d);

                if (magD > rootVSmall)
                {
                    d /= magD;

                    // Check how many empty directions are used by the edge.
                    label nEmptyDirs = 0;
                    label nNonEmptyDirs = 0;
                    for (direction cmpt=0; cmpt<vector::nComponents; cmpt++)
                    {
                        if (mag(d[cmpt]) > 1e-6)
                        {
                            if (directions[cmpt] == 0)
                            {
                                nEmptyDirs++;
                            }
                            else
                            {
                                nNonEmptyDirs++;
                            }
                        }
                    }

                    if (nEmptyDirs == 0)
                    {
                        // Purely in ok directions.
                    }
                    else if (nEmptyDirs == 1)
                    {
                        // Ok if purely in empty directions.
                        if (nNonEmptyDirs > 0)
                        {
                            edgesInError.insert(edge(p0, p1), facei);
                        }
                    }
                    else if (nEmptyDirs > 1)
                    {
                        // Always an error
                        edgesInError.insert(edge(p0, p1), facei);
                    }
                }
            }
        }
    }

    label nErrorEdges = returnReduce(edgesInError.size(), sumOp<label>());

    if (nErrorEdges > 0)
    {
        if (report)
        {
            Info<< " ***Number of edges not aligned with or perpendicular to "
                << "non-empty directions: " << nErrorEdges << endl;
        }

        if (setPtr)
        {
            setPtr->resize(2*edgesInError.size());
            forAllConstIter(EdgeMap<label>, edgesInError, iter)
            {
                setPtr->insert(iter.key()[0]);
                setPtr->insert(iter.key()[1]);
            }
        }

        return true;
    }
    else
    {
        if (report)
        {
            Info<< "    All edges aligned with or perpendicular to "
                << "non-empty directions." << endl;
        }
        return false;
    }
}


bool Foam::meshCheck::checkCellDeterminant
(
    const polyMesh& mesh,
    const bool report,
    labelHashSet* setPtr
)
{
    const vectorField& faceAreas = mesh.faceAreas();
    const Vector<label>& meshD = mesh.geometricD();

    const scalar warnDet = 1e-3;

    if (mesh.debug)
    {
        InfoInFunction << "Checking for under-determined cells" << endl;
    }

    tmp<scalarField> tcellDeterminant = meshCheck::cellDeterminant
    (
        mesh,
        meshD,
        faceAreas,
        syncTools::getInternalOrCoupledFaces(mesh)
    );
    scalarField& cellDeterminant = tcellDeterminant.ref();


    label nErrorCells = 0;
    scalar minDet = min(cellDeterminant);
    scalar sumDet = sum(cellDeterminant);

    forAll(cellDeterminant, celli)
    {
        if (cellDeterminant[celli] < warnDet)
        {
            if (setPtr)
            {
                setPtr->insert(celli);
            }

            nErrorCells++;
        }
    }

    reduce(nErrorCells, sumOp<label>());
    reduce(minDet, minOp<scalar>());
    reduce(sumDet, sumOp<scalar>());
    label nSummed = returnReduce(cellDeterminant.size(), sumOp<label>());

    if (report)
    {
        if (nSummed > 0)
        {
            Info<< "    Cell determinant (wellposedness) : minimum: " << minDet
                << " average: " << sumDet/nSummed
                << endl;
        }
    }

    if (nErrorCells > 0)
    {
        if (report)
        {
            Info<< " ***Cells with small determinant (< "
                << warnDet << ") found, number of cells: "
                << nErrorCells << endl;
        }

        return true;
    }
    else
    {
        if (report)
        {
            Info<< "    Cell determinant check OK." << endl;
        }

        return false;
    }

    return false;
}


bool Foam::meshCheck::checkFaceWeight
(
    const polyMesh& mesh,
    const bool report,
    const scalar minWeight,
    labelHashSet* setPtr
)
{
    if (mesh.debug)
    {
        InfoInFunction << "Checking for low face interpolation weights" << endl;
    }

    const vectorField& fCtrs = mesh.faceCentres();
    const vectorField& fAreas = mesh.faceAreas();
    const vectorField& cellCtrs = mesh.cellCentres();

    tmp<scalarField> tfaceWght = meshCheck::faceWeights
    (
        mesh,
        fCtrs,
        fAreas,
        cellCtrs
    );
    scalarField& faceWght = tfaceWght.ref();


    label nErrorFaces = 0;
    scalar minDet = great;
    scalar sumDet = 0.0;
    label nSummed = 0;

    // Statistics only for internal and masters of coupled faces
    PackedBoolList isMasterFace(syncTools::getInternalOrMasterFaces(mesh));

    forAll(faceWght, facei)
    {
        if (faceWght[facei] < minWeight)
        {
            // Note: insert both sides of coupled faces
            if (setPtr)
            {
                setPtr->insert(facei);
            }

            nErrorFaces++;
        }

        // Note: statistics only on master of coupled faces
        if (isMasterFace[facei])
        {
            minDet = min(minDet, faceWght[facei]);
            sumDet += faceWght[facei];
            nSummed++;
        }
    }

    reduce(nErrorFaces, sumOp<label>());
    reduce(minDet, minOp<scalar>());
    reduce(sumDet, sumOp<scalar>());
    reduce(nSummed, sumOp<label>());

    if (report)
    {
        if (nSummed > 0)
        {
            Info<< "    Face interpolation weight : minimum: " << minDet
                << " average: " << sumDet/nSummed
                << endl;
        }
    }

    if (nErrorFaces > 0)
    {
        if (report)
        {
            Info<< " ***Faces with small interpolation weight (< " << minWeight
                << ") found, number of faces: "
                << nErrorFaces << endl;
        }

        return true;
    }
    else
    {
        if (report)
        {
            Info<< "    Face interpolation weight check OK." << endl;
        }

        return false;
    }

    return false;
}


bool Foam::meshCheck::checkVolRatio
(
    const polyMesh& mesh,
    const bool report,
    const scalar minRatio,
    labelHashSet* setPtr
)
{
    if (mesh.debug)
    {
        InfoInFunction << "Checking for volume ratio < " << minRatio << endl;
    }

    const scalarField& cellVols = mesh.cellVolumes();

    tmp<scalarField> tvolRatio = meshCheck::volRatio(mesh, cellVols);
    scalarField& volRatio = tvolRatio.ref();


    label nErrorFaces = 0;
    scalar minDet = great;
    scalar sumDet = 0.0;
    label nSummed = 0;

    // Statistics only for internal and masters of coupled faces
    PackedBoolList isMasterFace(syncTools::getInternalOrMasterFaces(mesh));

    forAll(volRatio, facei)
    {
        if (volRatio[facei] < minRatio)
        {
            // Note: insert both sides of coupled faces
            if (setPtr)
            {
                setPtr->insert(facei);
            }

            nErrorFaces++;
        }

        // Note: statistics only on master of coupled faces
        if (isMasterFace[facei])
        {
            minDet = min(minDet, volRatio[facei]);
            sumDet += volRatio[facei];
            nSummed++;
        }
    }

    reduce(nErrorFaces, sumOp<label>());
    reduce(minDet, minOp<scalar>());
    reduce(sumDet, sumOp<scalar>());
    reduce(nSummed, sumOp<label>());

    if (report)
    {
        if (nSummed > 0)
        {
            Info<< "    Face volume ratio : minimum: " << minDet
                << " average: " << sumDet/nSummed
                << endl;
        }
    }

    if (nErrorFaces > 0)
    {
        if (report)
        {
            Info<< " ***Faces with small volume ratio (< " << minRatio
                << ") found, number of faces: "
                << nErrorFaces << endl;
        }

        return true;
    }
    else
    {
        if (report)
        {
            Info<< "    Face volume ratio check OK." << endl;
        }

        return false;
    }

    return false;
}


// ************************************************************************* //
