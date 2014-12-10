/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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

#include "primitiveMesh.H"
#include "pyramidPointFaceRef.H"
#include "ListOps.H"
#include "unitConversion.H"
#include "SortableList.H"
#include "EdgeMap.H"
#include "primitiveMeshTools.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

Foam::scalar Foam::primitiveMesh::closedThreshold_  = 1.0e-6;
Foam::scalar Foam::primitiveMesh::aspectThreshold_  = 1000;
Foam::scalar Foam::primitiveMesh::nonOrthThreshold_ = 70;    // deg
Foam::scalar Foam::primitiveMesh::skewThreshold_    = 4;
Foam::scalar Foam::primitiveMesh::planarCosAngle_   = 1.0e-6;


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::primitiveMesh::checkClosedBoundary
(
    const vectorField& areas,
    const bool report,
    const PackedBoolList& internalOrCoupledFaces
) const
{
    if (debug)
    {
        Info<< "bool primitiveMesh::checkClosedBoundary("
            << "const bool) const: "
            << "checking whether the boundary is closed" << endl;
    }

    // Loop through all boundary faces and sum up the face area vectors.
    // For a closed boundary, this should be zero in all vector components

    vector sumClosed(vector::zero);
    scalar sumMagClosedBoundary = 0;

    for (label faceI = nInternalFaces(); faceI < areas.size(); faceI++)
    {
        if (!internalOrCoupledFaces.size() || !internalOrCoupledFaces[faceI])
        {
            sumClosed += areas[faceI];
            sumMagClosedBoundary += mag(areas[faceI]);
        }
    }

    reduce(sumClosed, sumOp<vector>());
    reduce(sumMagClosedBoundary, sumOp<scalar>());

    vector openness = sumClosed/(sumMagClosedBoundary + VSMALL);

    if (cmptMax(cmptMag(openness)) > closedThreshold_)
    {
        if (debug || report)
        {
            Info<< " ***Boundary openness " << openness
                << " possible hole in boundary description."
                << endl;
        }

        return true;
    }
    else
    {
        if (debug || report)
        {
            Info<< "    Boundary openness " << openness << " OK."
                << endl;
        }

        return false;
    }
}


bool Foam::primitiveMesh::checkClosedCells
(
    const vectorField& faceAreas,
    const scalarField& cellVolumes,
    const bool report,
    labelHashSet* setPtr,
    labelHashSet* aspectSetPtr,
    const Vector<label>& meshD
) const
{
    if (debug)
    {
        Info<< "bool primitiveMesh::checkClosedCells("
            << "const bool, labelHashSet*, labelHashSet*"
            << ", const Vector<label>&) const: "
            << "checking whether cells are closed" << endl;
    }

    // Check that all cells labels are valid
    const cellList& c = cells();

    label nErrorClosed = 0;

    forAll(c, cI)
    {
        const cell& curCell = c[cI];

        if (min(curCell) < 0 || max(curCell) > nFaces())
        {
            if (setPtr)
            {
                setPtr->insert(cI);
            }

            nErrorClosed++;
        }
    }

    if (nErrorClosed > 0)
    {
        if (debug || report)
        {
            Info<< " ***Cells with invalid face labels found, number of cells "
                << nErrorClosed << endl;
        }

        return true;
    }


    scalarField openness;
    scalarField aspectRatio;
    primitiveMeshTools::cellClosedness
    (
        *this,
        meshD,
        faceAreas,
        cellVolumes,
        openness,
        aspectRatio
    );

    label nOpen = 0;
    scalar maxOpennessCell = max(openness);
    label nAspect = 0;
    scalar maxAspectRatio = max(aspectRatio);

    // Check the sums
    forAll(openness, cellI)
    {
        if (openness[cellI] > closedThreshold_)
        {
            if (setPtr)
            {
                setPtr->insert(cellI);
            }

            nOpen++;
        }

        if (aspectRatio[cellI] > aspectThreshold_)
        {
            if (aspectSetPtr)
            {
                aspectSetPtr->insert(cellI);
            }

            nAspect++;
        }
    }

    reduce(nOpen, sumOp<label>());
    reduce(maxOpennessCell, maxOp<scalar>());

    reduce(nAspect, sumOp<label>());
    reduce(maxAspectRatio, maxOp<scalar>());


    if (nOpen > 0)
    {
        if (debug || report)
        {
            Info<< " ***Open cells found, max cell openness: "
                << maxOpennessCell << ", number of open cells " << nOpen
                << endl;
        }

        return true;
    }

    if (nAspect > 0)
    {
        if (debug || report)
        {
            Info<< " ***High aspect ratio cells found, Max aspect ratio: "
                << maxAspectRatio
                << ", number of cells " << nAspect
                << endl;
        }

        return true;
    }

    if (debug || report)
    {
        Info<< "    Max cell openness = " << maxOpennessCell << " OK." << nl
            << "    Max aspect ratio = " << maxAspectRatio << " OK."
            << endl;
    }

    return false;
}


bool Foam::primitiveMesh::checkFaceAreas
(
    const vectorField& faceAreas,
    const bool report,
    const bool detailedReport,
    labelHashSet* setPtr
) const
{
    if (debug)
    {
        Info<< "bool primitiveMesh::checkFaceAreas("
            << "const bool, labelHashSet*) const: "
            << "checking face area magnitudes" << endl;
    }

    const scalarField magFaceAreas(mag(faceAreas));

    scalar minArea = GREAT;
    scalar maxArea = -GREAT;

    forAll(magFaceAreas, faceI)
    {
        if (magFaceAreas[faceI] < VSMALL)
        {
            if (setPtr)
            {
                setPtr->insert(faceI);
            }
            if (detailedReport)
            {
                if (isInternalFace(faceI))
                {
                    Pout<< "Zero or negative face area detected for "
                        << "internal face "<< faceI << " between cells "
                        << faceOwner()[faceI] << " and "
                        << faceNeighbour()[faceI]
                        << ".  Face area magnitude = " << magFaceAreas[faceI]
                        << endl;
                }
                else
                {
                    Pout<< "Zero or negative face area detected for "
                        << "boundary face " << faceI << " next to cell "
                        << faceOwner()[faceI] << ".  Face area magnitude = "
                        << magFaceAreas[faceI] << endl;
                }
            }
        }

        minArea = min(minArea, magFaceAreas[faceI]);
        maxArea = max(maxArea, magFaceAreas[faceI]);
    }

    reduce(minArea, minOp<scalar>());
    reduce(maxArea, maxOp<scalar>());

    if (minArea < VSMALL)
    {
        if (debug || report)
        {
            Info<< " ***Zero or negative face area detected.  "
                "Minimum area: " << minArea << endl;
        }

        return true;
    }
    else
    {
        if (debug || report)
        {
            Info<< "    Minimum face area = " << minArea
                << ". Maximum face area = " << maxArea
                << ".  Face area magnitudes OK." << endl;
        }

        return false;
    }
}


bool Foam::primitiveMesh::checkCellVolumes
(
    const scalarField& vols,
    const bool report,
    const bool detailedReport,
    labelHashSet* setPtr
) const
{
    if (debug)
    {
        Info<< "bool primitiveMesh::checkCellVolumes("
            << "const bool, labelHashSet*) const: "
            << "checking cell volumes" << endl;
    }

    scalar minVolume = GREAT;
    scalar maxVolume = -GREAT;

    label nNegVolCells = 0;

    forAll(vols, cellI)
    {
        if (vols[cellI] < VSMALL)
        {
            if (setPtr)
            {
                setPtr->insert(cellI);
            }
            if (detailedReport)
            {
                Pout<< "Zero or negative cell volume detected for cell "
                    << cellI << ".  Volume = " << vols[cellI] << endl;
            }

            nNegVolCells++;
        }

        minVolume = min(minVolume, vols[cellI]);
        maxVolume = max(maxVolume, vols[cellI]);
    }

    reduce(minVolume, minOp<scalar>());
    reduce(maxVolume, maxOp<scalar>());
    reduce(nNegVolCells, sumOp<label>());

    if (minVolume < VSMALL)
    {
        if (debug || report)
        {
            Info<< " ***Zero or negative cell volume detected.  "
                << "Minimum negative volume: " << minVolume
                << ", Number of negative volume cells: " << nNegVolCells
                << endl;
        }

        return true;
    }
    else
    {
        if (debug || report)
        {
            Info<< "    Min volume = " << minVolume
                << ". Max volume = " << maxVolume
                << ".  Total volume = " << gSum(vols)
                << ".  Cell volumes OK." << endl;
        }

        return false;
    }
}


bool Foam::primitiveMesh::checkFaceOrthogonality
(
    const vectorField& fAreas,
    const vectorField& cellCtrs,
    const bool report,
    labelHashSet* setPtr
) const
{
    if (debug)
    {
        Info<< "bool primitiveMesh::checkFaceOrthogonality("
            << "const bool, labelHashSet*) const: "
            << "checking mesh non-orthogonality" << endl;
    }


    tmp<scalarField> tortho = primitiveMeshTools::faceOrthogonality
    (
        *this,
        fAreas,
        cellCtrs
    );
    const scalarField& ortho = tortho();

    // Severe nonorthogonality threshold
    const scalar severeNonorthogonalityThreshold =
        ::cos(degToRad(nonOrthThreshold_));

    scalar minDDotS = min(ortho);

    scalar sumDDotS = sum(ortho);

    label severeNonOrth = 0;

    label errorNonOrth = 0;


    forAll(ortho, faceI)
    {
        if (ortho[faceI] < severeNonorthogonalityThreshold)
        {
            if (ortho[faceI] > SMALL)
            {
                if (setPtr)
                {
                    setPtr->insert(faceI);
                }

                severeNonOrth++;
            }
            else
            {
                if (setPtr)
                {
                    setPtr->insert(faceI);
                }

                errorNonOrth++;
            }
        }
    }

    reduce(minDDotS, minOp<scalar>());
    reduce(sumDDotS, sumOp<scalar>());
    reduce(severeNonOrth, sumOp<label>());
    reduce(errorNonOrth, sumOp<label>());

    if (debug || report)
    {
        label neiSize = ortho.size();
        reduce(neiSize, sumOp<label>());

        if (neiSize > 0)
        {
            if (debug || report)
            {
                Info<< "    Mesh non-orthogonality Max: "
                    << radToDeg(::acos(minDDotS))
                    << " average: " << radToDeg(::acos(sumDDotS/neiSize))
                    << endl;
            }
        }

        if (severeNonOrth > 0)
        {
            Info<< "   *Number of severely non-orthogonal faces: "
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


bool Foam::primitiveMesh::checkFacePyramids
(
    const pointField& points,
    const vectorField& ctrs,
    const bool report,
    const bool detailedReport,
    const scalar minPyrVol,
    labelHashSet* setPtr
) const
{
    if (debug)
    {
        Info<< "bool primitiveMesh::checkFacePyramids("
            << "const bool, const scalar, labelHashSet*) const: "
            << "checking face orientation" << endl;
    }

    const labelList& own = faceOwner();
    const labelList& nei = faceNeighbour();
    const faceList& f = faces();


    scalarField ownPyrVol;
    scalarField neiPyrVol;
    primitiveMeshTools::facePyramidVolume
    (
        *this,
        points,
        ctrs,
        ownPyrVol,
        neiPyrVol
    );


    label nErrorPyrs = 0;

    forAll(ownPyrVol, faceI)
    {
        if (ownPyrVol[faceI] < minPyrVol)
        {
            if (setPtr)
            {
                setPtr->insert(faceI);
            }
            if (detailedReport)
            {
                Pout<< "Negative pyramid volume: " << ownPyrVol[faceI]
                    << " for face " << faceI << " " << f[faceI]
                    << "  and owner cell: " << own[faceI] << endl
                    << "Owner cell vertex labels: "
                    << cells()[own[faceI]].labels(faces())
                    << endl;
            }

            nErrorPyrs++;
        }

        if (isInternalFace(faceI))
        {
            if (neiPyrVol[faceI] < minPyrVol)
            {
                if (setPtr)
                {
                    setPtr->insert(faceI);
                }
                if (detailedReport)
                {
                    Pout<< "Negative pyramid volume: " << neiPyrVol[faceI]
                        << " for face " << faceI << " " << f[faceI]
                        << "  and neighbour cell: " << nei[faceI] << nl
                        << "Neighbour cell vertex labels: "
                        << cells()[nei[faceI]].labels(faces())
                        << endl;
                }
                nErrorPyrs++;
            }
        }
    }

    reduce(nErrorPyrs, sumOp<label>());

    if (nErrorPyrs > 0)
    {
        if (debug || report)
        {
            Info<< " ***Error in face pyramids: "
                << nErrorPyrs << " faces are incorrectly oriented."
                << endl;
        }

        return true;
    }
    else
    {
        if (debug || report)
        {
            Info<< "    Face pyramids OK." << endl;
        }

        return false;
    }
}


bool Foam::primitiveMesh::checkFaceSkewness
(
    const pointField& points,
    const vectorField& fCtrs,
    const vectorField& fAreas,
    const vectorField& cellCtrs,
    const bool report,
    labelHashSet* setPtr
) const
{
    if (debug)
    {
        Info<< "bool primitiveMesh::checkFaceSkewnesss("
            << "const bool, labelHashSet*) const: "
            << "checking face skewness" << endl;
    }

    // Warn if the skew correction vector is more than skewWarning times
    // larger than the face area vector

    tmp<scalarField> tskewness = primitiveMeshTools::faceSkewness
    (
        *this,
        points,
        fCtrs,
        fAreas,
        cellCtrs
    );
    const scalarField& skewness = tskewness();

    scalar maxSkew = max(skewness);
    label nWarnSkew = 0;

    forAll(skewness, faceI)
    {
        // Check if the skewness vector is greater than the PN vector.
        // This does not cause trouble but is a good indication of a poor mesh.
        if (skewness[faceI] > skewThreshold_)
        {
            if (setPtr)
            {
                setPtr->insert(faceI);
            }

            nWarnSkew++;
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


// Check convexity of angles in a face. Allow a slight non-convexity.
// E.g. maxDeg = 10 allows for angles < 190 (or 10 degrees concavity)
// (if truly concave and points not visible from face centre the face-pyramid
//  check in checkMesh will fail)
bool Foam::primitiveMesh::checkFaceAngles
(
    const pointField& points,
    const vectorField& faceAreas,
    const bool report,
    const scalar maxDeg,
    labelHashSet* setPtr
) const
{
    if (debug)
    {
        Info<< "bool primitiveMesh::checkFaceAngles"
            << "(const bool, const scalar, labelHashSet*) const: "
            << "checking face angles" << endl;
    }

    if (maxDeg < -SMALL || maxDeg > 180+SMALL)
    {
        FatalErrorIn
        (
            "primitiveMesh::checkFaceAngles"
            "(const bool, const scalar, labelHashSet*)"
        )   << "maxDeg should be [0..180] but is now " << maxDeg
            << exit(FatalError);
    }

    const scalar maxSin = Foam::sin(degToRad(maxDeg));


    tmp<scalarField> tfaceAngles = primitiveMeshTools::faceConcavity
    (
        maxSin,
        *this,
        points,
        faceAreas
    );
    const scalarField& faceAngles = tfaceAngles();

    scalar maxEdgeSin = max(faceAngles);

    label nConcave = 0;

    forAll(faceAngles, faceI)
    {
        if (faceAngles[faceI] > SMALL)
        {
            nConcave++;

            if (setPtr)
            {
                setPtr->insert(faceI);
            }
        }
    }

    reduce(nConcave, sumOp<label>());
    reduce(maxEdgeSin, maxOp<scalar>());

    if (nConcave > 0)
    {
        scalar maxConcaveDegr =
            radToDeg(Foam::asin(Foam::min(1.0, maxEdgeSin)));

        if (debug || report)
        {
            Info<< "   *There are " << nConcave
                << " faces with concave angles between consecutive"
                << " edges. Max concave angle = " << maxConcaveDegr
                << " degrees." << endl;
        }

        return true;
    }
    else
    {
        if (debug || report)
        {
            Info<< "    All angles in faces OK." << endl;
        }

        return false;
    }
}


bool Foam::primitiveMesh::checkFaceFlatness
(
    const pointField& points,
    const vectorField& faceCentres,
    const vectorField& faceAreas,
    const bool report,
    const scalar warnFlatness,
    labelHashSet* setPtr
) const
{
    if (debug)
    {
        Info<< "bool primitiveMesh::checkFaceFlatness"
            << "(const bool, const scalar, labelHashSet*) const: "
            << "checking face flatness" << endl;
    }

    if (warnFlatness < 0 || warnFlatness > 1)
    {
        FatalErrorIn
        (
            "primitiveMesh::checkFaceFlatness"
            "(const bool, const scalar, labelHashSet*)"
        )   << "warnFlatness should be [0..1] but is now " << warnFlatness
            << exit(FatalError);
    }

    const faceList& fcs = faces();

    tmp<scalarField> tfaceFlatness = primitiveMeshTools::faceFlatness
    (
        *this,
        points,
        faceCentres,
        faceAreas
    );
    const scalarField& faceFlatness = tfaceFlatness();

    scalarField magAreas(mag(faceAreas));

    scalar minFlatness = GREAT;
    scalar sumFlatness = 0;
    label nSummed = 0;
    label nWarped = 0;

    forAll(faceFlatness, faceI)
    {
        if (fcs[faceI].size() > 3 && magAreas[faceI] > VSMALL)
        {
            sumFlatness += faceFlatness[faceI];
            nSummed++;

            minFlatness = min(minFlatness, faceFlatness[faceI]);

            if (faceFlatness[faceI] < warnFlatness)
            {
                nWarped++;

                if (setPtr)
                {
                    setPtr->insert(faceI);
                }
            }
        }
    }


    reduce(nWarped, sumOp<label>());
    reduce(minFlatness, minOp<scalar>());

    reduce(nSummed, sumOp<label>());
    reduce(sumFlatness, sumOp<scalar>());

    if (debug || report)
    {
        if (nSummed > 0)
        {
            Info<< "    Face flatness (1 = flat, 0 = butterfly) : min = "
                << minFlatness << "  average = " << sumFlatness / nSummed
                << endl;
        }
    }


    if (nWarped> 0)
    {
        if (debug || report)
        {
            Info<< "   *There are " << nWarped
                << " faces with ratio between projected and actual area < "
                << warnFlatness << endl;

            Info<< "    Minimum ratio (minimum flatness, maximum warpage) = "
                << minFlatness << endl;
        }

        return true;
    }
    else
    {
        if (debug || report)
        {
            Info<< "    All face flatness OK." << endl;
        }

        return false;
    }
}


bool Foam::primitiveMesh::checkConcaveCells
(
    const vectorField& fAreas,
    const pointField& fCentres,
    const bool report,
    labelHashSet* setPtr
) const
{
    if (debug)
    {
        Info<< "bool primitiveMesh::checkConcaveCells(const bool"
            << ", labelHashSet*) const: "
            << "checking for concave cells" << endl;
    }

    const cellList& c = cells();
    const labelList& fOwner = faceOwner();

    label nConcaveCells = 0;

    forAll(c, cellI)
    {
        const cell& cFaces = c[cellI];

        bool concave = false;

        forAll(cFaces, i)
        {
            if (concave)
            {
                break;
            }

            label fI = cFaces[i];

            const point& fC = fCentres[fI];

            vector fN = fAreas[fI];

            fN /= max(mag(fN), VSMALL);

            // Flip normal if required so that it is always pointing out of
            // the cell
            if (fOwner[fI] != cellI)
            {
                fN *= -1;
            }

            // Is the centre of any other face of the cell on the
            // wrong side of the plane of this face?

            forAll(cFaces, j)
            {
                if (j != i)
                {
                    label fJ = cFaces[j];

                    const point& pt = fCentres[fJ];

                    // If the cell is concave, the point will be on the
                    // positive normal side of the plane of f, defined by
                    // its centre and normal, and the angle between (pt -
                    // fC) and fN will be less than 90 degrees, so the dot
                    // product will be positive.

                    vector pC = (pt - fC);

                    pC /= max(mag(pC), VSMALL);

                    if ((pC & fN) > -planarCosAngle_)
                    {
                        // Concave or planar face

                        concave = true;

                        if (setPtr)
                        {
                            setPtr->insert(cellI);
                        }

                        nConcaveCells++;

                        break;
                    }
                }
            }
        }
    }

    reduce(nConcaveCells, sumOp<label>());

    if (nConcaveCells > 0)
    {
        if (debug || report)
        {
            Info<< " ***Concave cells (using face planes) found,"
                << " number of cells: " << nConcaveCells << endl;
        }

        return true;
    }
    else
    {
        if (debug || report)
        {
            Info<< "    Concave cell check OK." << endl;
        }

        return false;
    }

    return false;
}


// Topological tests

bool Foam::primitiveMesh::checkUpperTriangular
(
    const bool report,
    labelHashSet* setPtr
) const
{
    if (debug)
    {
        Info<< "bool primitiveMesh::checkUpperTriangular("
            << "const bool, labelHashSet*) const: "
            << "checking face ordering" << endl;
    }

    // Check whether internal faces are ordered in the upper triangular order
    const labelList& own = faceOwner();
    const labelList& nei = faceNeighbour();

    const cellList& c = cells();

    label internal = nInternalFaces();

    // Has error occurred?
    bool error = false;
    // Have multiple faces been detected?
    label nMultipleCells = false;

    // Loop through faceCells once more and make sure that for internal cell
    // the first label is smaller
    for (label faceI = 0; faceI < internal; faceI++)
    {
        if (own[faceI] >= nei[faceI])
        {
            error  = true;

            if (setPtr)
            {
                setPtr->insert(faceI);
            }
        }
    }

    // Loop through all cells. For each cell, find the face that is internal
    // and add it to the check list (upper triangular order).
    // Once the list is completed, check it against the faceCell list

    forAll(c, cellI)
    {
        const labelList& curFaces = c[cellI];

        // Neighbouring cells
        SortableList<label> nbr(curFaces.size());

        forAll(curFaces, i)
        {
            label faceI = curFaces[i];

            if (faceI >= nInternalFaces())
            {
                // Sort last
                nbr[i] = labelMax;
            }
            else
            {
                label nbrCellI = nei[faceI];

                if (nbrCellI == cellI)
                {
                    nbrCellI = own[faceI];
                }

                if (cellI < nbrCellI)
                {
                    // cellI is master
                    nbr[i] = nbrCellI;
                }
                else
                {
                    // nbrCell is master. Let it handle this face.
                    nbr[i] = labelMax;
                }
            }
        }

        nbr.sort();

        // Now nbr holds the cellCells in incremental order. Check:
        // - neighbouring cells appear only once. Since nbr is sorted this
        //   is simple check on consecutive elements
        // - faces indexed in same order as nbr are incrementing as well.

        label prevCell = nbr[0];
        label prevFace = curFaces[nbr.indices()[0]];

        bool hasMultipleFaces = false;

        for (label i = 1; i < nbr.size(); i++)
        {
            label thisCell = nbr[i];
            label thisFace = curFaces[nbr.indices()[i]];

            if (thisCell == labelMax)
            {
                break;
            }

            if (thisCell == prevCell)
            {
                hasMultipleFaces = true;

                if (setPtr)
                {
                    setPtr->insert(prevFace);
                    setPtr->insert(thisFace);
                }
            }
            else if (thisFace < prevFace)
            {
                error = true;

                if (setPtr)
                {
                    setPtr->insert(thisFace);
                }
            }

            prevCell = thisCell;
            prevFace = thisFace;
        }

        if (hasMultipleFaces)
        {
            nMultipleCells++;
        }
    }

    reduce(error, orOp<bool>());
    reduce(nMultipleCells, sumOp<label>());

    if ((debug || report) && nMultipleCells > 0)
    {
        Info<< "  <<Found " << nMultipleCells
            << " neighbouring cells with multiple inbetween faces." << endl;
    }

    if (error)
    {
        if (debug || report)
        {
            Info<< " ***Faces not in upper triangular order." << endl;
        }

        return true;
    }
    else
    {
        if (debug || report)
        {
            Info<< "    Upper triangular ordering OK." << endl;
        }

        return false;
    }
}


bool Foam::primitiveMesh::checkCellsZipUp
(
    const bool report,
    labelHashSet* setPtr
) const
{
    if (debug)
    {
        Info<< "bool primitiveMesh::checkCellsZipUp("
            << "const bool, labelHashSet*) const: "
            << "checking topological cell openness" << endl;
    }

    label nOpenCells = 0;

    const faceList& f = faces();
    const cellList& c = cells();

    forAll(c, cellI)
    {
        const labelList& curFaces = c[cellI];

        const edgeList cellEdges = c[cellI].edges(f);

        labelList edgeUsage(cellEdges.size(), 0);

        forAll(curFaces, faceI)
        {
            edgeList curFaceEdges = f[curFaces[faceI]].edges();

            forAll(curFaceEdges, faceEdgeI)
            {
                const edge& curEdge = curFaceEdges[faceEdgeI];

                forAll(cellEdges, cellEdgeI)
                {
                    if (cellEdges[cellEdgeI] == curEdge)
                    {
                        edgeUsage[cellEdgeI]++;
                        break;
                    }
                }
            }
        }

        edgeList singleEdges(cellEdges.size());
        label nSingleEdges = 0;

        forAll(edgeUsage, edgeI)
        {
            if (edgeUsage[edgeI] == 1)
            {
                singleEdges[nSingleEdges] = cellEdges[edgeI];
                nSingleEdges++;
            }
            else if (edgeUsage[edgeI] != 2)
            {
                if (setPtr)
                {
                    setPtr->insert(cellI);
                }
            }
        }

        if (nSingleEdges > 0)
        {
            if (setPtr)
            {
                setPtr->insert(cellI);
            }

            nOpenCells++;
        }
    }

    reduce(nOpenCells, sumOp<label>());

    if (nOpenCells > 0)
    {
        if (debug || report)
        {
            Info<< " ***Open cells found, number of cells: " << nOpenCells
                << ". This problem may be fixable using the zipUpMesh utility."
                << endl;
        }

        return true;
    }
    else
    {
        if (debug || report)
        {
            Info<< "    Topological cell zip-up check OK." << endl;
        }

        return false;
    }
}


// Vertices of face within point range and unique.
bool Foam::primitiveMesh::checkFaceVertices
(
    const bool report,
    labelHashSet* setPtr
) const
{
    if (debug)
    {
        Info<< "bool primitiveMesh::checkFaceVertices("
            << "const bool, labelHashSet*) const: "
            << "checking face vertices" << endl;
    }

    // Check that all vertex labels are valid
    const faceList& f = faces();

    label nErrorFaces = 0;

    forAll(f, fI)
    {
        const face& curFace = f[fI];

        if (min(curFace) < 0 || max(curFace) > nPoints())
        {
            if (setPtr)
            {
                setPtr->insert(fI);
            }

            nErrorFaces++;
        }

        // Uniqueness of vertices
        labelHashSet facePoints(2*curFace.size());

        forAll(curFace, fp)
        {
            bool inserted = facePoints.insert(curFace[fp]);

            if (!inserted)
            {
                if (setPtr)
                {
                    setPtr->insert(fI);
                }

                nErrorFaces++;
            }
        }
    }

    reduce(nErrorFaces, sumOp<label>());

    if (nErrorFaces > 0)
    {
        if (debug || report)
        {
            Info<< "    Faces with invalid vertex labels found, "
                << " number of faces: " << nErrorFaces << endl;
        }

        return true;
    }
    else
    {
        if (debug || report)
        {
            Info<< "    Face vertices OK." << endl;
        }

        return false;
    }
}


bool Foam::primitiveMesh::checkPoints
(
    const bool report,
    labelHashSet* setPtr
) const
{
    if (debug)
    {
        Info<< "bool primitiveMesh::checkPoints"
            << "(const bool, labelHashSet*) const: "
            << "checking points" << endl;
    }

    label nFaceErrors = 0;
    label nCellErrors = 0;

    const labelListList& pf = pointFaces();

    forAll(pf, pointI)
    {
        if (pf[pointI].empty())
        {
            if (setPtr)
            {
                setPtr->insert(pointI);
            }

            nFaceErrors++;
        }
    }


    forAll(pf, pointI)
    {
        const labelList& pc = pointCells(pointI);

        if (pc.empty())
        {
            if (setPtr)
            {
                setPtr->insert(pointI);
            }

            nCellErrors++;
        }
    }

    reduce(nFaceErrors, sumOp<label>());
    reduce(nCellErrors, sumOp<label>());

    if (nFaceErrors > 0 || nCellErrors > 0)
    {
        if (debug || report)
        {
            Info<< " ***Unused points found in the mesh, "
                   "number unused by faces: " << nFaceErrors
                << " number unused by cells: " << nCellErrors
                << endl;
        }

        return true;
    }
    else
    {
        if (debug || report)
        {
            Info<< "    Point usage OK." << endl;
        }

        return false;
    }
}


// Check if all points on face are shared between faces.
bool Foam::primitiveMesh::checkDuplicateFaces
(
    const label faceI,
    const Map<label>& nCommonPoints,
    label& nBaffleFaces,
    labelHashSet* setPtr
) const
{
    bool error = false;

    forAllConstIter(Map<label>, nCommonPoints, iter)
    {
        label nbFaceI = iter.key();
        label nCommon = iter();

        const face& curFace = faces()[faceI];
        const face& nbFace = faces()[nbFaceI];

        if (nCommon == nbFace.size() || nCommon == curFace.size())
        {
            if (nbFace.size() != curFace.size())
            {
                error = true;
            }
            else
            {
                nBaffleFaces++;
            }

            if (setPtr)
            {
                setPtr->insert(faceI);
                setPtr->insert(nbFaceI);
            }
        }
    }

    return error;
}


// Check that shared points are in consecutive order.
bool Foam::primitiveMesh::checkCommonOrder
(
    const label faceI,
    const Map<label>& nCommonPoints,
    labelHashSet* setPtr
) const
{
    bool error = false;

    forAllConstIter(Map<label>, nCommonPoints, iter)
    {
        label nbFaceI = iter.key();
        label nCommon = iter();

        const face& curFace = faces()[faceI];
        const face& nbFace = faces()[nbFaceI];

        if
        (
            nCommon >= 2
         && nCommon != nbFace.size()
         && nCommon != curFace.size()
        )
        {
            forAll(curFace, fp)
            {
                // Get the index in the neighbouring face shared with curFace
                label nb = findIndex(nbFace, curFace[fp]);

                if (nb != -1)
                {

                    // Check the whole face from nb onwards for shared vertices
                    // with neighbouring face. Rule is that any shared vertices
                    // should be consecutive on both faces i.e. if they are
                    // vertices fp,fp+1,fp+2 on one face they should be
                    // vertices nb, nb+1, nb+2 (or nb+2, nb+1, nb) on the
                    // other face.


                    // Vertices before and after on curFace
                    label fpPlus1 = curFace.fcIndex(fp);
                    label fpMin1  = curFace.rcIndex(fp);

                    // Vertices before and after on nbFace
                    label nbPlus1 = nbFace.fcIndex(nb);
                    label nbMin1  = nbFace.rcIndex(nb);

                    // Find order of walking by comparing next points on both
                    // faces.
                    label curInc = labelMax;
                    label nbInc = labelMax;

                    if (nbFace[nbPlus1] == curFace[fpPlus1])
                    {
                        curInc = 1;
                        nbInc = 1;
                    }
                    else if (nbFace[nbPlus1] == curFace[fpMin1])
                    {
                        curInc = -1;
                        nbInc = 1;
                    }
                    else if (nbFace[nbMin1] == curFace[fpMin1])
                    {
                        curInc = -1;
                        nbInc = -1;
                    }
                    else
                    {
                        curInc = 1;
                        nbInc = -1;
                    }


                    // Pass1: loop until start of common vertices found.
                    label curNb = nb;
                    label curFp = fp;

                    do
                    {
                        curFp += curInc;

                        if (curFp >= curFace.size())
                        {
                            curFp = 0;
                        }
                        else if (curFp < 0)
                        {
                            curFp = curFace.size()-1;
                        }

                        curNb += nbInc;

                        if (curNb >= nbFace.size())
                        {
                            curNb = 0;
                        }
                        else if (curNb < 0)
                        {
                            curNb = nbFace.size()-1;
                        }
                    } while (curFace[curFp] == nbFace[curNb]);


                    // Pass2: check equality walking from curFp, curNb
                    // in opposite order.

                    curInc = -curInc;
                    nbInc = -nbInc;

                    for (label commonI = 0; commonI < nCommon; commonI++)
                    {
                        curFp += curInc;

                        if (curFp >= curFace.size())
                        {
                            curFp = 0;
                        }
                        else if (curFp < 0)
                        {
                            curFp = curFace.size()-1;
                        }

                        curNb += nbInc;

                        if (curNb >= nbFace.size())
                        {
                            curNb = 0;
                        }
                        else if (curNb < 0)
                        {
                            curNb = nbFace.size()-1;
                        }

                        if (curFace[curFp] != nbFace[curNb])
                        {
                            if (setPtr)
                            {
                                setPtr->insert(faceI);
                                setPtr->insert(nbFaceI);
                            }

                            error = true;

                            break;
                        }
                    }


                    // Done the curFace - nbFace combination.
                    break;
                }
            }
        }
    }

    return error;
}


// Checks common vertices between faces. If more than 2 they should be
// consecutive on both faces.
bool Foam::primitiveMesh::checkFaceFaces
(
    const bool report,
    labelHashSet* setPtr
) const
{
    if (debug)
    {
        Info<< "bool primitiveMesh::checkFaceFaces(const bool, labelHashSet*)"
            << " const: " << "checking face-face connectivity" << endl;
    }

    const labelListList& pf = pointFaces();

    label nBaffleFaces = 0;
    label nErrorDuplicate = 0;
    label nErrorOrder = 0;
    Map<label> nCommonPoints(100);

    for (label faceI = 0; faceI < nFaces(); faceI++)
    {
        const face& curFace = faces()[faceI];

        // Calculate number of common points between current faceI and
        // neighbouring face. Store on map.
        nCommonPoints.clear();

        forAll(curFace, fp)
        {
            label pointI = curFace[fp];

            const labelList& nbs = pf[pointI];

            forAll(nbs, nbI)
            {
                label nbFaceI = nbs[nbI];

                if (faceI < nbFaceI)
                {
                    // Only check once for each combination of two faces.

                    Map<label>::iterator fnd = nCommonPoints.find(nbFaceI);

                    if (fnd == nCommonPoints.end())
                    {
                        // First common vertex found.
                        nCommonPoints.insert(nbFaceI, 1);
                    }
                    else
                    {
                        fnd()++;
                    }
                }
            }
        }

        // Perform various checks on common points

        // Check all vertices shared (duplicate point)
        if (checkDuplicateFaces(faceI, nCommonPoints, nBaffleFaces, setPtr))
        {
            nErrorDuplicate++;
        }

        // Check common vertices are consecutive on both faces
        if (checkCommonOrder(faceI, nCommonPoints, setPtr))
        {
            nErrorOrder++;
        }
    }

    reduce(nBaffleFaces, sumOp<label>());
    reduce(nErrorDuplicate, sumOp<label>());
    reduce(nErrorOrder, sumOp<label>());

    if (nBaffleFaces)
    {
        Info<< "    Number of identical duplicate faces (baffle faces): "
            << nBaffleFaces << endl;
    }

    if (nErrorDuplicate > 0 || nErrorOrder > 0)
    {
        // These are actually warnings, not errors.
        if (nErrorDuplicate > 0)
        {
            Info<< "  <<Number of duplicate (not baffle) faces found: "
                << nErrorDuplicate
                << ". This might indicate a problem." << endl;
        }

        if (nErrorOrder > 0)
        {
            Info<< "  <<Number of faces with non-consecutive shared points: "
                << nErrorOrder << ". This might indicate a problem." << endl;
        }

        return false;   //return true;
    }
    else
    {
        if (debug || report)
        {
            Info<< "    Face-face connectivity OK." << endl;
        }

        return false;
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::primitiveMesh::checkClosedBoundary(const bool report) const
{
    return checkClosedBoundary(faceAreas(), report, PackedBoolList(0));
}


bool Foam::primitiveMesh::checkClosedCells
(
    const bool report,
    labelHashSet* setPtr,
    labelHashSet* aspectSetPtr,
     const Vector<label>& solutionD
) const
{
    return checkClosedCells
    (
        faceAreas(),
        cellVolumes(),
        report,
        setPtr,
        aspectSetPtr,
        solutionD
    );
}


bool Foam::primitiveMesh::checkFaceAreas
(
    const bool report,
    labelHashSet* setPtr
) const
{
    return checkFaceAreas
    (
        faceAreas(),
        report,
        false,  // detailedReport,
        setPtr
    );
}


bool Foam::primitiveMesh::checkCellVolumes
(
    const bool report,
    labelHashSet* setPtr
) const
{
    return checkCellVolumes
    (
        cellVolumes(),
        report,
        false,  // detailedReport,
        setPtr
    );
}


bool Foam::primitiveMesh::checkFaceOrthogonality
(
    const bool report,
    labelHashSet* setPtr
) const
{
    return checkFaceOrthogonality
    (
        faceAreas(),
        cellCentres(),
        report,
        setPtr
    );
}


bool Foam::primitiveMesh::checkFacePyramids
(
    const bool report,
    const scalar minPyrVol,
    labelHashSet* setPtr
) const
{
    return checkFacePyramids
    (
        points(),
        cellCentres(),
        report,
        false,  // detailedReport,
        minPyrVol,
        setPtr
    );
}


bool Foam::primitiveMesh::checkFaceSkewness
(
    const bool report,
    labelHashSet* setPtr
) const
{
    return checkFaceSkewness
    (
        points(),
        faceCentres(),
        faceAreas(),
        cellCentres(),
        report,
        setPtr
    );
}


bool Foam::primitiveMesh::checkFaceAngles
(
    const bool report,
    const scalar maxDeg,
    labelHashSet* setPtr
) const
{
    return checkFaceAngles
    (
        points(),
        faceAreas(),
        report,
        maxDeg,
        setPtr
    );
}


bool Foam::primitiveMesh::checkFaceFlatness
(
    const bool report,
    const scalar warnFlatness,
    labelHashSet* setPtr
) const
{
    return checkFaceFlatness
    (
        points(),
        faceCentres(),
        faceAreas(),
        report,
        warnFlatness,
        setPtr
    );
}


bool Foam::primitiveMesh::checkConcaveCells
(
    const bool report,
    labelHashSet* setPtr
) const
{
    return checkConcaveCells
    (
        faceAreas(),
        faceCentres(),
        report,
        setPtr
    );
}


bool Foam::primitiveMesh::checkTopology(const bool report) const
{
    label noFailedChecks = 0;

    if (checkPoints(report)) noFailedChecks++;
    if (checkUpperTriangular(report)) noFailedChecks++;
    if (checkCellsZipUp(report)) noFailedChecks++;
    if (checkFaceVertices(report)) noFailedChecks++;
    if (checkFaceFaces(report)) noFailedChecks++;

    if (noFailedChecks == 0)
    {
        if (debug || report)
        {
            Info<< "    Mesh topology OK." << endl;
        }

        return false;
    }
    else
    {
        if (debug || report)
        {
            Info<< "    Failed " << noFailedChecks
                << " mesh topology checks." << endl;
        }

        return true;
    }
}


bool Foam::primitiveMesh::checkGeometry(const bool report) const
{
    label noFailedChecks = 0;

    if (checkClosedBoundary(report)) noFailedChecks++;
    if (checkClosedCells(report)) noFailedChecks++;
    if (checkFaceAreas(report)) noFailedChecks++;
    if (checkCellVolumes(report)) noFailedChecks++;
    if (checkFaceOrthogonality(report)) noFailedChecks++;
    if (checkFacePyramids(report)) noFailedChecks++;
    if (checkFaceSkewness(report)) noFailedChecks++;

    if (noFailedChecks == 0)
    {
        if (debug || report)
        {
            Info<< "    Mesh geometry OK." << endl;
        }
        return false;
    }
    else
    {
        if (debug || report)
        {
            Info<< "    Failed " << noFailedChecks
                << " mesh geometry checks." << endl;
        }

        return true;
    }
}


bool Foam::primitiveMesh::checkMesh(const bool report) const
{
    if (debug)
    {
        Info<< "bool primitiveMesh::checkMesh(const bool report) const: "
            << "checking primitiveMesh" << endl;
    }

    label noFailedChecks = checkTopology(report) + checkGeometry(report);

    if (noFailedChecks == 0)
    {
        if (debug || report)
        {
            Info<< "Mesh OK." << endl;
        }

        return false;
    }
    else
    {
        if (debug || report)
        {
            Info<< "    Failed " << noFailedChecks
                << " mesh checks." << endl;
        }

        return true;
    }
}


Foam::scalar Foam::primitiveMesh::setClosedThreshold(const scalar val)
{
    scalar prev = closedThreshold_;
    closedThreshold_ = val;

    return prev;
}


Foam::scalar Foam::primitiveMesh::setAspectThreshold(const scalar val)
{
    scalar prev = aspectThreshold_;
    aspectThreshold_ = val;

    return prev;
}


Foam::scalar Foam::primitiveMesh::setNonOrthThreshold(const scalar val)
{
    scalar prev = nonOrthThreshold_;
    nonOrthThreshold_ = val;

    return prev;
}


Foam::scalar Foam::primitiveMesh::setSkewThreshold(const scalar val)
{
    scalar prev = skewThreshold_;
    skewThreshold_ = val;

    return prev;
}


// ************************************************************************* //
