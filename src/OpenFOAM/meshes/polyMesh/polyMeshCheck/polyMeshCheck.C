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

#include "polyMeshCheck.H"
#include "polyMeshTools.H"
#include "unitConversion.H"
#include "syncTools.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

Foam::scalar Foam::polyMeshCheck::closedThreshold  = 1.0e-6;
Foam::scalar Foam::polyMeshCheck::aspectThreshold  = 1000;
Foam::scalar Foam::polyMeshCheck::nonOrthThreshold = 70;    // deg
Foam::scalar Foam::polyMeshCheck::skewThreshold    = 4;
Foam::scalar Foam::polyMeshCheck::planarCosAngle   = 1.0e-6;

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::polyMesh::checkFaceOrthogonality
(
    const bool report,
    labelHashSet* setPtr
) const
{
    if (debug)
    {
        InfoInFunction << "Checking mesh non-orthogonality" << endl;
    }

    const vectorField& fAreas = faceAreas();
    const vectorField& cellCtrs = cellCentres();

    // Calculate orthogonality for all internal and coupled boundary faces
    // (1 for uncoupled boundary faces)
    tmp<scalarField> tortho = polyMeshTools::faceOrthogonality
    (
        *this,
        fAreas,
        cellCtrs
    );
    const scalarField& ortho = tortho.ref();

    // Severe nonorthogonality threshold
    const scalar severeNonorthogonalityThreshold =
        ::cos(degToRad(polyMeshCheck::nonOrthThreshold));


    scalar minDDotS = great;
    scalar sumDDotS = 0.0;
    label nSummed = 0;
    label severeNonOrth = 0;
    label errorNonOrth = 0;


    // Statistics only for internal and masters of coupled faces
    PackedBoolList isMasterFace(syncTools::getInternalOrMasterFaces(*this));

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

    if (debug || report)
    {
        if (nSummed > 0)
        {
            if (debug || report)
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
                << polyMeshCheck::nonOrthThreshold << " degrees) faces: "
                << severeNonOrth << "." << endl;
        }
    }

    if (errorNonOrth > 0)
    {
        if (debug || report)
        {
            Info<< " ***Number of non-orthogonality errors: "
                << errorNonOrth << "." << endl;
        }

        return true;
    }
    else
    {
        if (debug || report)
        {
            Info<< "    Non-orthogonality check OK." << endl;
        }

        return false;
    }
}


bool Foam::polyMesh::checkFaceSkewness
(
    const bool report,
    labelHashSet* setPtr
) const
{
    if (debug)
    {
        InfoInFunction << "Checking face skewness" << endl;
    }

    const pointField& points = this->points();
    const vectorField& fCtrs = faceCentres();
    const vectorField& fAreas = faceAreas();
    const vectorField& cellCtrs = cellCentres();

    // Warn if the skew correction vector is more than skewWarning times
    // larger than the face area vector

    tmp<scalarField> tskew = polyMeshTools::faceSkewness
    (
        *this,
        points,
        fCtrs,
        fAreas,
        cellCtrs
    );
    const scalarField& skew = tskew.ref();

    scalar maxSkew = max(skew);
    label nWarnSkew = 0;

    // Statistics only for all faces except slave coupled faces
    PackedBoolList isMasterFace(syncTools::getMasterFaces(*this));

    forAll(skew, facei)
    {
        // Check if the skewness vector is greater than the PN vector.
        // This does not cause trouble but is a good indication of a poor mesh.
        if (skew[facei] > polyMeshCheck::skewThreshold)
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
        if (debug || report)
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
        if (debug || report)
        {
            Info<< "    Max skewness = " << maxSkew << " OK." << endl;
        }

        return false;
    }
}


bool Foam::polyMesh::checkEdgeAlignment
(
    const bool report,
    const Vector<label>& directions,
    labelHashSet* setPtr
) const
{
    // Check 1D/2Dness of edges. Gets passed the non-empty directions and
    // checks all edges in the mesh whether they:
    // - have no component in a non-empty direction or
    // - are only in a single non-empty direction.
    // Empty direction info is passed in as a vector of labels (synchronised)
    // which are 1 if the direction is non-empty, 0 if it is.

    if (debug)
    {
        InfoInFunction << "Checking edge alignment" << endl;
    }

    const pointField& p = points();

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


    const faceList& fcs = faces();

    EdgeMap<label> edgesInError;

    forAll(fcs, facei)
    {
        const face& f = fcs[facei];

        forAll(f, fp)
        {
            label p0 = f[fp];
            label p1 = f.nextLabel(fp);
            if (p0 < p1)
            {
                vector d(p[p1]-p[p0]);
                scalar magD = mag(d);

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
        if (debug || report)
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
        if (debug || report)
        {
            Info<< "    All edges aligned with or perpendicular to "
                << "non-empty directions." << endl;
        }
        return false;
    }
}


bool Foam::polyMesh::checkCellDeterminant
(
    const bool report,
    labelHashSet* setPtr
) const
{
    const vectorField& faceAreas = this->faceAreas();
    const Vector<label>& meshD = geometricD();

    const scalar warnDet = 1e-3;

    if (debug)
    {
        InfoInFunction << "Checking for under-determined cells" << endl;
    }

    tmp<scalarField> tcellDeterminant = primitiveMeshTools::cellDeterminant
    (
        *this,
        meshD,
        faceAreas,
        syncTools::getInternalOrCoupledFaces(*this)
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

    if (debug || report)
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
        if (debug || report)
        {
            Info<< " ***Cells with small determinant (< "
                << warnDet << ") found, number of cells: "
                << nErrorCells << endl;
        }

        return true;
    }
    else
    {
        if (debug || report)
        {
            Info<< "    Cell determinant check OK." << endl;
        }

        return false;
    }

    return false;
}


bool Foam::polyMesh::checkFaceWeight
(
    const bool report,
    const scalar minWeight,
    labelHashSet* setPtr
) const
{
    if (debug)
    {
        InfoInFunction << "Checking for low face interpolation weights" << endl;
    }

    const vectorField& fCtrs = faceCentres();
    const vectorField& fAreas = faceAreas();
    const vectorField& cellCtrs = this->cellCentres();

    tmp<scalarField> tfaceWght = polyMeshTools::faceWeights
    (
        *this,
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
    PackedBoolList isMasterFace(syncTools::getInternalOrMasterFaces(*this));

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

    if (debug || report)
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
        if (debug || report)
        {
            Info<< " ***Faces with small interpolation weight (< " << minWeight
                << ") found, number of faces: "
                << nErrorFaces << endl;
        }

        return true;
    }
    else
    {
        if (debug || report)
        {
            Info<< "    Face interpolation weight check OK." << endl;
        }

        return false;
    }

    return false;
}


bool Foam::polyMesh::checkVolRatio
(
    const bool report,
    const scalar minRatio,
    labelHashSet* setPtr
) const
{
    if (debug)
    {
        InfoInFunction << "Checking for volume ratio < " << minRatio << endl;
    }

    const scalarField& cellVols = cellVolumes();

    tmp<scalarField> tvolRatio = polyMeshTools::volRatio(*this, cellVols);
    scalarField& volRatio = tvolRatio.ref();


    label nErrorFaces = 0;
    scalar minDet = great;
    scalar sumDet = 0.0;
    label nSummed = 0;

    // Statistics only for internal and masters of coupled faces
    PackedBoolList isMasterFace(syncTools::getInternalOrMasterFaces(*this));

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

    if (debug || report)
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
        if (debug || report)
        {
            Info<< " ***Faces with small volume ratio (< " << minRatio
                << ") found, number of faces: "
                << nErrorFaces << endl;
        }

        return true;
    }
    else
    {
        if (debug || report)
        {
            Info<< "    Face volume ratio check OK." << endl;
        }

        return false;
    }

    return false;
}



bool Foam::polyMeshCheck::checkTopology(const polyMesh& mesh, const bool report)
{
    label noFailedChecks = 0;

    if (mesh.checkPoints(report)) noFailedChecks++;
    if (mesh.checkUpperTriangular(report)) noFailedChecks++;
    if (mesh.checkCellsZipUp(report)) noFailedChecks++;
    if (mesh.checkFaceVertices(report)) noFailedChecks++;
    if (mesh.checkFaceFaces(report)) noFailedChecks++;

    if (noFailedChecks == 0)
    {
        if (report)
        {
            Info<< "    Mesh topology OK." << endl;
        }

        return false;
    }
    else
    {
        if (report)
        {
            Info<< "    Failed " << noFailedChecks
                << " mesh topology checks." << endl;
        }

        return true;
    }
}


bool Foam::polyMeshCheck::checkGeometry(const polyMesh& mesh, const bool report)
{
    label noFailedChecks = 0;

    if (mesh.checkClosedBoundary(report)) noFailedChecks++;
    if (mesh.checkClosedCells(report)) noFailedChecks++;
    if (mesh.checkFaceAreas(report)) noFailedChecks++;
    if (mesh.checkCellVolumes(report)) noFailedChecks++;
    if (mesh.checkFaceOrthogonality(report)) noFailedChecks++;
    if (mesh.checkFacePyramids(report)) noFailedChecks++;
    if (mesh.checkFaceSkewness(report)) noFailedChecks++;

    if (noFailedChecks == 0)
    {
        if (report)
        {
            Info<< "    Mesh geometry OK." << endl;
        }

        return false;
    }
    else
    {
        if (report)
        {
            Info<< "    Failed " << noFailedChecks
                << " mesh geometry checks." << endl;
        }

        return true;
    }
}


bool Foam::polyMeshCheck::checkMesh(const polyMesh& mesh, const bool report)
{
    const label noFailedChecks =
        checkTopology(mesh, report)
      + checkGeometry(mesh, report);

    if (noFailedChecks == 0)
    {
        if (report)
        {
            Info<< "Mesh OK." << endl;
        }

        return false;
    }
    else
    {
        if (report)
        {
            Info<< "    Failed " << noFailedChecks
                << " mesh checks." << endl;
        }

        return true;
    }
}


// ************************************************************************* //
