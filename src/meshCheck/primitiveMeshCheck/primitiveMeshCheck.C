/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2024 OpenFOAM Foundation
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
#include "pyramidPointFaceRef.H"
#include "PackedBoolList.H"
#include "unitConversion.H"
#include "SortableList.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::meshCheck::faceSkewness
(
    const primitiveMesh& mesh,
    const pointField& p,
    const vectorField& fCtrs,
    const vectorField& fAreas,

    const label facei,
    const point& ownCc,
    const point& neiCc
)
{
    const vector Cpf = fCtrs[facei] - ownCc;
    const vector d = neiCc - ownCc;

    // Skewness vector
    const vector sv =
        Cpf
      - ((fAreas[facei] & Cpf)/((fAreas[facei] & d) + rootVSmall))*d;
    const vector svHat = sv/(mag(sv) + rootVSmall);

    // Normalisation distance calculated as the approximate distance
    // from the face centre to the edge of the face in the direction
    // of the skewness
    scalar fd = 0.2*mag(d) + rootVSmall;
    const face& f = mesh.faces()[facei];
    forAll(f, pi)
    {
        fd = max(fd, mag(svHat & (p[f[pi]] - fCtrs[facei])));
    }

    // Normalised skewness
    return mag(sv)/fd;
}


Foam::scalar Foam::meshCheck::boundaryFaceSkewness
(
    const primitiveMesh& mesh,
    const pointField& p,
    const vectorField& fCtrs,
    const vectorField& fAreas,

    const label facei,
    const point& ownCc
)
{
    const vector Cpf = fCtrs[facei] - ownCc;

    vector normal = fAreas[facei];
    normal /= mag(normal) + rootVSmall;
    const vector d = normal*(normal & Cpf);


    // Skewness vector
    const vector sv =
        Cpf
      - ((fAreas[facei] & Cpf)/((fAreas[facei] & d) + rootVSmall))*d;
    const vector svHat = sv/(mag(sv) + rootVSmall);

    // Normalisation distance calculated as the approximate distance
    // from the face centre to the edge of the face in the direction
    // of the skewness
    scalar fd = 0.4*mag(d) + rootVSmall;
    const face& f = mesh.faces()[facei];
    forAll(f, pi)
    {
        fd = max(fd, mag(svHat & (p[f[pi]] - fCtrs[facei])));
    }

    // Normalised skewness
    return mag(sv)/fd;
}


Foam::scalar Foam::meshCheck::faceOrthogonality
(
    const point& ownCc,
    const point& neiCc,
    const vector& s
)
{
    const vector d = neiCc - ownCc;
    return (d & s)/(mag(d)*mag(s) + rootVSmall);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::scalarField> Foam::meshCheck::faceOrthogonality
(
    const primitiveMesh& mesh,
    const vectorField& areas,
    const vectorField& cc
)
{
    const labelList& own = mesh.faceOwner();
    const labelList& nei = mesh.faceNeighbour();

    tmp<scalarField> tortho(new scalarField(mesh.nInternalFaces()));
    scalarField& ortho = tortho.ref();

    // Internal faces
    forAll(nei, facei)
    {
        ortho[facei] = faceOrthogonality
        (
            cc[own[facei]],
            cc[nei[facei]],
            areas[facei]
        );
    }

    return tortho;
}


Foam::tmp<Foam::scalarField> Foam::meshCheck::faceSkewness
(
    const primitiveMesh& mesh,
    const pointField& p,
    const vectorField& fCtrs,
    const vectorField& fAreas,
    const vectorField& cellCtrs
)
{
    const labelList& own = mesh.faceOwner();
    const labelList& nei = mesh.faceNeighbour();

    tmp<scalarField> tskew(new scalarField(mesh.nFaces()));
    scalarField& skew = tskew.ref();

    forAll(nei, facei)
    {
        skew[facei] = faceSkewness
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

    for (label facei = mesh.nInternalFaces(); facei < mesh.nFaces(); facei++)
    {
        skew[facei] = boundaryFaceSkewness
        (
            mesh,
            p,
            fCtrs,
            fAreas,
            facei,
            cellCtrs[own[facei]]
        );
    }

    return tskew;
}


void Foam::meshCheck::facePyramidVolume
(
    const primitiveMesh& mesh,
    const pointField& points,
    const vectorField& ctrs,

    scalarField& ownPyrVol,
    scalarField& neiPyrVol
)
{
    const labelList& own = mesh.faceOwner();
    const labelList& nei = mesh.faceNeighbour();
    const faceList& f = mesh.faces();

    ownPyrVol.setSize(mesh.nFaces());
    neiPyrVol.setSize(mesh.nInternalFaces());

    forAll(f, facei)
    {
        // Create the owner pyramid
        ownPyrVol[facei] = -pyramidPointFaceRef
        (
            f[facei],
            ctrs[own[facei]]
        ).mag(points);

        if (mesh.isInternalFace(facei))
        {
            // Create the neighbour pyramid - it will have positive volume
            neiPyrVol[facei] = pyramidPointFaceRef
            (
                f[facei],
                ctrs[nei[facei]]
            ).mag(points);
        }
    }
}


void Foam::meshCheck::cellClosedness
(
    const primitiveMesh& mesh,
    const Vector<label>& meshD,
    const vectorField& areas,
    const scalarField& vols,

    scalarField& openness,
    scalarField& aratio
)
{
    const labelList& own = mesh.faceOwner();
    const labelList& nei = mesh.faceNeighbour();

    // Loop through cell faces and sum up the face area vectors for each cell.
    // This should be zero in all vector components

    vectorField sumClosed(mesh.nCells(), Zero);
    vectorField sumMagClosed(mesh.nCells(), Zero);

    forAll(own, facei)
    {
        // Add to owner
        sumClosed[own[facei]] += areas[facei];
        sumMagClosed[own[facei]] += cmptMag(areas[facei]);
    }

    forAll(nei, facei)
    {
        // Subtract from neighbour
        sumClosed[nei[facei]] -= areas[facei];
        sumMagClosed[nei[facei]] += cmptMag(areas[facei]);
    }


    label nDims = 0;
    for (direction dir = 0; dir < vector::nComponents; dir++)
    {
        if (meshD[dir] == 1)
        {
            nDims++;
        }
    }


    // Check the sums
    openness.setSize(mesh.nCells());
    aratio.setSize(mesh.nCells());

    forAll(sumClosed, celli)
    {
        scalar maxOpenness = 0;

        for (direction cmpt=0; cmpt<vector::nComponents; cmpt++)
        {
            maxOpenness = max
            (
                maxOpenness,
                mag(sumClosed[celli][cmpt])
               /(sumMagClosed[celli][cmpt] + rootVSmall)
            );
        }
        openness[celli] = maxOpenness;

        // Calculate the aspect ration as the maximum of Cartesian component
        // aspect ratio to the total area hydraulic area aspect ratio
        scalar minCmpt = vGreat;
        scalar maxCmpt = -vGreat;
        for (direction dir = 0; dir < vector::nComponents; dir++)
        {
            if (meshD[dir] == 1)
            {
                minCmpt = min(minCmpt, sumMagClosed[celli][dir]);
                maxCmpt = max(maxCmpt, sumMagClosed[celli][dir]);
            }
        }

        scalar aspectRatio = maxCmpt/(minCmpt + rootVSmall);
        if (nDims == 3)
        {
            scalar v = max(rootVSmall, vols[celli]);

            aspectRatio = max
            (
                aspectRatio,
                1.0/6.0*cmptSum(sumMagClosed[celli])/pow(v, 2.0/3.0)
            );
        }

        aratio[celli] = aspectRatio;
    }
}


Foam::tmp<Foam::scalarField> Foam::meshCheck::faceConcavity
(
    const scalar maxSin,
    const primitiveMesh& mesh,
    const pointField& p,
    const vectorField& faceAreas
)
{
    const faceList& fcs = mesh.faces();

    vectorField faceNormals(faceAreas);
    faceNormals /= mag(faceNormals) + rootVSmall;

    tmp<scalarField> tfaceAngles(new scalarField(mesh.nFaces()));
    scalarField& faceAngles = tfaceAngles.ref();


    forAll(fcs, facei)
    {
        const face& f = fcs[facei];

        // Get edge from f[0] to f[size-1];
        vector ePrev(p[f.first()] - p[f.last()]);
        scalar magEPrev = mag(ePrev);
        ePrev /= magEPrev + rootVSmall;

        scalar maxEdgeSin = 0.0;

        forAll(f, fp0)
        {
            // Get vertex after fp
            const label fp1 = f.fcIndex(fp0);

            // Normalised vector between two consecutive points
            vector e10(p[f[fp1]] - p[f[fp0]]);
            const scalar magE10 = mag(e10);
            e10 /= magE10 + rootVSmall;

            if (magEPrev > small && magE10 > small)
            {
                vector edgeNormal = ePrev ^ e10;
                const scalar magEdgeNormal = mag(edgeNormal);

                if (magEdgeNormal < maxSin)
                {
                    // Edges (almost) aligned -> face is ok.
                }
                else
                {
                    // Check normal
                    edgeNormal /= magEdgeNormal;

                    if ((edgeNormal & faceNormals[facei]) < small)
                    {
                        maxEdgeSin = max(maxEdgeSin, magEdgeNormal);
                    }
                }
            }

            ePrev = e10;
            magEPrev = magE10;
        }

        faceAngles[facei] = maxEdgeSin;
    }

    return tfaceAngles;
}


Foam::tmp<Foam::scalarField> Foam::meshCheck::faceFlatness
(
    const primitiveMesh& mesh,
    const pointField& p,
    const vectorField& fCtrs,
    const vectorField& faceAreas
)
{
    const faceList& fcs = mesh.faces();

    // Areas are calculated as the sum of areas. (see
    // primitiveMeshFaceCentresAndAreas.C)
    scalarField magAreas(mag(faceAreas));

    tmp<scalarField> tfaceFlatness(new scalarField(mesh.nFaces(), 1.0));
    scalarField& faceFlatness = tfaceFlatness.ref();


    forAll(fcs, facei)
    {
        const face& f = fcs[facei];

        if (f.size() > 3 && magAreas[facei] > rootVSmall)
        {
            const point& fc = fCtrs[facei];

            // Calculate the sum of magnitude of areas and compare to magnitude
            // of sum of areas.

            scalar sumA = 0.0;

            forAll(f, fp)
            {
                const point& thisPoint = p[f[fp]];
                const point& nextPoint = p[f.nextLabel(fp)];

                // Triangle around fc.
                const vector n = 0.5*((nextPoint - thisPoint)^(fc - thisPoint));
                sumA += mag(n);
            }

            faceFlatness[facei] = magAreas[facei]/(sumA + rootVSmall);
        }
    }

    return tfaceFlatness;
}


Foam::tmp<Foam::scalarField> Foam::meshCheck::cellDeterminant
(
    const primitiveMesh& mesh,
    const Vector<label>& meshD,
    const vectorField& faceAreas,
    const PackedBoolList& internalOrCoupledFace
)
{
    // Determine number of dimensions and (for 2D) missing dimension
    label nDims = 0;
    label twoD = -1;
    for (direction dir = 0; dir < vector::nComponents; dir++)
    {
        if (meshD[dir] == 1)
        {
            nDims++;
        }
        else
        {
            twoD = dir;
        }
    }

    tmp<scalarField> tcellDeterminant(new scalarField(mesh.nCells()));
    scalarField& cellDeterminant = tcellDeterminant.ref();

    const cellList& c = mesh.cells();

    if (nDims == 1)
    {
        cellDeterminant = 1.0;
    }
    else
    {
        forAll(c, celli)
        {
            const labelList& curFaces = c[celli];

            // Calculate local normalisation factor
            scalar avgArea = 0;

            label nInternalFaces = 0;

            forAll(curFaces, i)
            {
                if (internalOrCoupledFace[curFaces[i]])
                {
                    avgArea += mag(faceAreas[curFaces[i]]);

                    nInternalFaces++;
                }
            }

            if (nInternalFaces == 0)
            {
                cellDeterminant[celli] = 0;
            }
            else
            {
                avgArea /= nInternalFaces;

                symmTensor areaTensor(Zero);

                forAll(curFaces, i)
                {
                    if (internalOrCoupledFace[curFaces[i]])
                    {
                        areaTensor += sqr(faceAreas[curFaces[i]]/avgArea);
                    }
                }

                if (nDims == 2)
                {
                    // Add the missing eigenvector (such that it does not
                    // affect the determinant)
                    if (twoD == 0)
                    {
                        areaTensor.xx() = 1;
                    }
                    else if (twoD == 1)
                    {
                        areaTensor.yy() = 1;
                    }
                    else
                    {
                        areaTensor.zz() = 1;
                    }
                }

                cellDeterminant[celli] = mag(det(areaTensor));
            }
        }
    }

    return tcellDeterminant;
}


bool Foam::meshCheck::checkClosedBoundary
(
    const primitiveMesh& mesh,
    const scalar closedThreshold,
    const bool report
)
{
    if (mesh.debug)
    {
        InfoInFunction
            << "Checking whether the boundary is closed" << endl;
    }

    const vectorField& areas = mesh.faceAreas();

    // Loop through all boundary faces and sum up the face area vectors.
    // For a closed boundary, this should be zero in all vector components

    vector sumClosed(Zero);
    scalar sumMagClosedBoundary = 0;

    for (label facei = mesh.nInternalFaces(); facei < areas.size(); facei++)
    {
        sumClosed += areas[facei];
        sumMagClosedBoundary += mag(areas[facei]);
    }

    reduce(sumClosed, sumOp<vector>());
    reduce(sumMagClosedBoundary, sumOp<scalar>());

    vector openness = sumClosed/(sumMagClosedBoundary + vSmall);

    if (cmptMax(cmptMag(openness)) > closedThreshold)
    {
        if (report)
        {
            Info<< " ***Boundary openness " << openness
                << " possible hole in boundary description."
                << endl;
        }

        return true;
    }
    else
    {
        if (report)
        {
            Info<< "    Boundary openness " << openness << " OK."
                << endl;
        }

        return false;
    }
}


bool Foam::meshCheck::checkClosedCells
(
    const primitiveMesh& mesh,
    const scalar closedThreshold,
    const scalar aspectThreshold,
    const bool report,
    labelHashSet* setPtr,
    labelHashSet* aspectSetPtr,
    const Vector<label>& meshD
)
{
    if (mesh.debug)
    {
        InfoInFunction
            << "Checking whether cells are closed" << endl;
    }

    const vectorField& faceAreas = mesh.faceAreas();
    const scalarField& cellVolumes = mesh.cellVolumes();

    // Check that all cells labels are valid
    const cellList& c = mesh.cells();

    label nErrorClosed = 0;

    forAll(c, cI)
    {
        const cell& curCell = c[cI];

        if (min(curCell) < 0 || max(curCell) > mesh.nFaces())
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
        if (report)
        {
            Info<< " ***Cells with invalid face labels found, number of cells "
                << nErrorClosed << endl;
        }

        return true;
    }


    scalarField openness;
    scalarField aspectRatio;
    meshCheck::cellClosedness
    (
        mesh,
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
    forAll(openness, celli)
    {
        if (openness[celli] > closedThreshold)
        {
            if (setPtr)
            {
                setPtr->insert(celli);
            }

            nOpen++;
        }

        if (aspectRatio[celli] > aspectThreshold)
        {
            if (aspectSetPtr)
            {
                aspectSetPtr->insert(celli);
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
        if (report)
        {
            Info<< " ***Open cells found, max cell openness: "
                << maxOpennessCell << ", number of open cells " << nOpen
                << endl;
        }

        return true;
    }

    if (nAspect > 0)
    {
        if (report)
        {
            Info<< " ***High aspect ratio cells found, Max aspect ratio: "
                << maxAspectRatio
                << ", number of cells " << nAspect
                << endl;
        }

        return true;
    }

    if (report)
    {
        Info<< "    Max cell openness = " << maxOpennessCell << " OK." << nl
            << "    Max aspect ratio = " << maxAspectRatio << " OK."
            << endl;
    }

    return false;
}


bool Foam::meshCheck::checkFaceAreas
(
    const primitiveMesh& mesh,
    const bool report,
    labelHashSet* setPtr
)
{
    if (mesh.debug)
    {
        InfoInFunction << "Checking face area magnitudes" << endl;
    }

    const vectorField& faceAreas = mesh.faceAreas();
    const scalarField magFaceAreas(mag(faceAreas));

    scalar minArea = great;
    scalar maxArea = -great;

    forAll(magFaceAreas, facei)
    {
        if (magFaceAreas[facei] < vSmall)
        {
            if (setPtr)
            {
                setPtr->insert(facei);
            }
        }

        minArea = min(minArea, magFaceAreas[facei]);
        maxArea = max(maxArea, magFaceAreas[facei]);
    }

    reduce(minArea, minOp<scalar>());
    reduce(maxArea, maxOp<scalar>());

    if (minArea < vSmall)
    {
        if (report)
        {
            Info<< " ***Zero or negative face area detected.  "
                "Minimum area: " << minArea << endl;
        }

        return true;
    }
    else
    {
        if (report)
        {
            Info<< "    Minimum face area = " << minArea
                << ". Maximum face area = " << maxArea
                << ".  Face area magnitudes OK." << endl;
        }

        return false;
    }
}


bool Foam::meshCheck::checkCellVolumes
(
    const primitiveMesh& mesh,
    const bool report,
    labelHashSet* setPtr
)
{
    if (mesh.debug)
    {
        InfoInFunction << "Checking cell volumes" << endl;
    }

    const scalarField& vols = mesh.cellVolumes();
    scalar minVolume = great;
    scalar maxVolume = -great;

    label nNegVolCells = 0;

    forAll(vols, celli)
    {
        if (vols[celli] < vSmall)
        {
            if (setPtr)
            {
                setPtr->insert(celli);
            }

            nNegVolCells++;
        }

        minVolume = min(minVolume, vols[celli]);
        maxVolume = max(maxVolume, vols[celli]);
    }

    reduce(minVolume, minOp<scalar>());
    reduce(maxVolume, maxOp<scalar>());
    reduce(nNegVolCells, sumOp<label>());

    if (minVolume < vSmall)
    {
        if (report)
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
        if (report)
        {
            Info<< "    Min volume = " << minVolume
                << ". Max volume = " << maxVolume
                << ".  Total volume = " << gSum(vols)
                << ".  Cell volumes OK." << endl;
        }

        return false;
    }
}


bool Foam::meshCheck::checkFacePyramids
(
    const primitiveMesh& mesh,
    const bool report,
    const scalar minPyrVol,
    labelHashSet* setPtr
)
{
    if (mesh.debug)
    {
        InfoInFunction << "Checking face orientation" << endl;
    }

    const pointField& points = mesh.points();
    const vectorField& ctrs = mesh.cellCentres();


    scalarField ownPyrVol;
    scalarField neiPyrVol;
    meshCheck::facePyramidVolume
    (
        mesh,
        points,
        ctrs,
        ownPyrVol,
        neiPyrVol
    );


    label nErrorPyrs = 0;

    forAll(ownPyrVol, facei)
    {
        if (ownPyrVol[facei] < minPyrVol)
        {
            if (setPtr)
            {
                setPtr->insert(facei);
            }

            nErrorPyrs++;
        }

        if (mesh.isInternalFace(facei))
        {
            if (neiPyrVol[facei] < minPyrVol)
            {
                if (setPtr)
                {
                    setPtr->insert(facei);
                }
                nErrorPyrs++;
            }
        }
    }

    reduce(nErrorPyrs, sumOp<label>());

    if (nErrorPyrs > 0)
    {
        if (report)
        {
            Info<< " ***Error in face pyramids: "
                << nErrorPyrs << " faces are incorrectly oriented."
                << endl;
        }

        return true;
    }
    else
    {
        if (report)
        {
            Info<< "    Face pyramids OK." << endl;
        }

        return false;
    }
}


bool Foam::meshCheck::checkFaceAngles
(
    const primitiveMesh& mesh,
    const bool report,
    const scalar maxConcave,
    labelHashSet* setPtr
)
{
    if (mesh.debug)
    {
        InfoInFunction << "Checking face angles" << endl;
    }

    if (maxConcave < -small || maxConcave > degToRad(180)+small)
    {
        FatalErrorInFunction
            << "maxConcave should be [0..180] degrees but is "
            << radToDeg(maxConcave) << abort(FatalError);
    }

    const scalar maxSin = Foam::sin(maxConcave);

    const pointField& points = mesh.points();
    const vectorField& faceAreas = mesh.faceAreas();

    tmp<scalarField> tfaceAngles = meshCheck::faceConcavity
    (
        maxSin,
        mesh,
        points,
        faceAreas
    );
    const scalarField& faceAngles = tfaceAngles();

    scalar maxEdgeSin = max(faceAngles);

    label nConcave = 0;

    forAll(faceAngles, facei)
    {
        if (faceAngles[facei] > small)
        {
            nConcave++;

            if (setPtr)
            {
                setPtr->insert(facei);
            }
        }
    }

    reduce(nConcave, sumOp<label>());
    reduce(maxEdgeSin, maxOp<scalar>());

    if (nConcave > 0)
    {
        if (report)
        {
            Info<< "   *There are " << nConcave
                << " faces with concave angles between consecutive"
                << " edges. Max concave angle = "
                << radToDeg(Foam::asin(Foam::min(1.0, maxEdgeSin)))
                << " degrees." << endl;
        }

        return true;
    }
    else
    {
        if (report)
        {
            Info<< "    All angles in faces OK." << endl;
        }

        return false;
    }
}


bool Foam::meshCheck::checkFaceFlatness
(
    const primitiveMesh& mesh,
    const bool report,
    const scalar warnFlatness,
    labelHashSet* setPtr
)
{
    if (mesh.debug)
    {
        InfoInFunction << "Checking face flatness" << endl;
    }

    if (warnFlatness < 0 || warnFlatness > 1)
    {
        FatalErrorInFunction
            << "warnFlatness should be [0..1] but is now " << warnFlatness
            << exit(FatalError);
    }

    const pointField& points = mesh.points();
    const vectorField& faceCentres = mesh.faceCentres();
    const vectorField& faceAreas = mesh.faceAreas();
    const faceList& fcs = mesh.faces();

    tmp<scalarField> tfaceFlatness = meshCheck::faceFlatness
    (
        mesh,
        points,
        faceCentres,
        faceAreas
    );
    const scalarField& faceFlatness = tfaceFlatness();

    scalarField magAreas(mag(faceAreas));

    scalar minFlatness = great;
    scalar sumFlatness = 0;
    label nSummed = 0;
    label nWarped = 0;

    forAll(faceFlatness, facei)
    {
        if (fcs[facei].size() > 3 && magAreas[facei] > vSmall)
        {
            sumFlatness += faceFlatness[facei];
            nSummed++;

            minFlatness = min(minFlatness, faceFlatness[facei]);

            if (faceFlatness[facei] < warnFlatness)
            {
                nWarped++;

                if (setPtr)
                {
                    setPtr->insert(facei);
                }
            }
        }
    }


    reduce(nWarped, sumOp<label>());
    reduce(minFlatness, minOp<scalar>());

    reduce(nSummed, sumOp<label>());
    reduce(sumFlatness, sumOp<scalar>());

    if (report)
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
        if (report)
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
        if (report)
        {
            Info<< "    All face flatness OK." << endl;
        }

        return false;
    }
}


bool Foam::meshCheck::checkConcaveCells
(
    const primitiveMesh& mesh,
    const scalar planarCosAngle,
    const bool report,
    labelHashSet* setPtr
)
{
    if (mesh.debug)
    {
        InfoInFunction << "Checking for concave cells" << endl;
    }

    const vectorField& fAreas = mesh.faceAreas();
    const pointField& fCentres = mesh.faceCentres();
    const cellList& c = mesh.cells();
    const labelList& fOwner = mesh.faceOwner();

    label nConcaveCells = 0;

    forAll(c, celli)
    {
        const cell& cFaces = c[celli];

        bool concave = false;

        forAll(cFaces, i)
        {
            if (concave)
            {
                break;
            }

            const label fI = cFaces[i];

            const point& fC = fCentres[fI];

            vector fN = fAreas[fI];
            fN /= max(mag(fN), vSmall);

            // Flip normal if required so that it is always pointing out of
            // the cell
            if (fOwner[fI] != celli)
            {
                fN *= -1;
            }

            // Is the centre of any other face of the cell on the
            // wrong side of the plane of this face?

            forAll(cFaces, j)
            {
                if (j != i)
                {
                    const label fJ = cFaces[j];

                    const point& pt = fCentres[fJ];

                    // If the cell is concave, the point will be on the
                    // positive normal side of the plane of f, defined by
                    // its centre and normal, and the angle between (pt -
                    // fC) and fN will be less than 90 degrees, so the dot
                    // product will be positive.

                    vector pC = (pt - fC);
                    pC /= max(mag(pC), vSmall);

                    if ((pC & fN) > -planarCosAngle)
                    {
                        // Concave or planar face

                        concave = true;

                        if (setPtr)
                        {
                            setPtr->insert(celli);
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
        if (report)
        {
            Info<< " ***Concave cells (using face planes) found,"
                << " number of cells: " << nConcaveCells << endl;
        }

        return true;
    }
    else
    {
        if (report)
        {
            Info<< "    Concave cell check OK." << endl;
        }

        return false;
    }

    return false;
}


bool Foam::meshCheck::checkUpperTriangular
(
    const primitiveMesh& mesh,
    const bool report,
    labelHashSet* setPtr
)
{
    if (mesh.debug)
    {
        InfoInFunction << "Checking face ordering" << endl;
    }

    // Check whether internal faces are ordered in the upper triangular order
    const labelList& own = mesh.faceOwner();
    const labelList& nei = mesh.faceNeighbour();
    const cellList& c = mesh.cells();
    const label internal = mesh.nInternalFaces();

    // Has error occurred?
    bool error = false;

    // Have multiple faces been detected?
    label nMultipleCells = false;

    // Loop through faceCells once more and make sure that for internal cell
    // the first label is smaller
    for (label facei = 0; facei < internal; facei++)
    {
        if (own[facei] >= nei[facei])
        {
            error  = true;

            if (setPtr)
            {
                setPtr->insert(facei);
            }
        }
    }

    // Loop through all cells. For each cell, find the face that is internal
    // and add it to the check list (upper triangular order).
    // Once the list is completed, check it against the faceCell list

    forAll(c, celli)
    {
        const labelList& curFaces = c[celli];

        // Neighbouring cells
        SortableList<label> nbr(curFaces.size());

        forAll(curFaces, i)
        {
            label facei = curFaces[i];

            if (facei >= mesh.nInternalFaces())
            {
                // Sort last
                nbr[i] = labelMax;
            }
            else
            {
                label nbrCelli = nei[facei];

                if (nbrCelli == celli)
                {
                    nbrCelli = own[facei];
                }

                if (celli < nbrCelli)
                {
                    // celli is master
                    nbr[i] = nbrCelli;
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
            const label thisCell = nbr[i];
            const label thisFace = curFaces[nbr.indices()[i]];

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

    if ((report) && nMultipleCells > 0)
    {
        Info<< "  <<Found " << nMultipleCells
            << " neighbouring cells with multiple in between faces." << endl;
    }

    if (error)
    {
        if (report)
        {
            Info<< " ***Faces not in upper triangular order." << endl;
        }

        return true;
    }
    else
    {
        if (report)
        {
            Info<< "    Upper triangular ordering OK." << endl;
        }

        return false;
    }
}


bool Foam::meshCheck::checkCellsZipUp
(
    const primitiveMesh& mesh,
    const bool report,
    labelHashSet* setPtr
)
{
    if (mesh.debug)
    {
        InfoInFunction << "Checking topological cell openness" << endl;
    }

    label nOpenCells = 0;

    const faceList& f = mesh.faces();
    const cellList& c = mesh.cells();

    forAll(c, celli)
    {
        const labelList& curFaces = c[celli];

        const edgeList cellEdges = c[celli].edges(f);

        labelList edgeUsage(cellEdges.size(), 0);

        forAll(curFaces, facei)
        {
            const edgeList curFaceEdges = f[curFaces[facei]].edges();

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
                    setPtr->insert(celli);
                }
            }
        }

        if (nSingleEdges > 0)
        {
            if (setPtr)
            {
                setPtr->insert(celli);
            }

            nOpenCells++;
        }
    }

    reduce(nOpenCells, sumOp<label>());

    if (nOpenCells > 0)
    {
        if (report)
        {
            Info<< " ***Open cells found, number of cells: " << nOpenCells
                << ". This problem may be fixable using the zipUpMesh utility."
                << endl;
        }

        return true;
    }
    else
    {
        if (report)
        {
            Info<< "    Topological cell zip-up check OK." << endl;
        }

        return false;
    }
}


bool Foam::meshCheck::checkFaceVertices
(
    const primitiveMesh& mesh,
    const bool report,
    labelHashSet* setPtr
)
{
    if (mesh.debug)
    {
        InfoInFunction << "Checking face vertices" << endl;
    }

    // Check that all vertex labels are valid
    const faceList& f = mesh.faces();

    label nErrorFaces = 0;

    forAll(f, fI)
    {
        const face& curFace = f[fI];

        if (min(curFace) < 0 || max(curFace) > mesh.nPoints())
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
            const bool inserted = facePoints.insert(curFace[fp]);

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
        if (report)
        {
            Info<< " ***Faces with invalid vertex labels found, "
                << " number of faces: " << nErrorFaces << endl;
        }

        return true;
    }
    else
    {
        if (report)
        {
            Info<< "    Face vertices OK." << endl;
        }

        return false;
    }
}


bool Foam::meshCheck::checkPoints
(
    const primitiveMesh& mesh,
    const bool report,
    labelHashSet* setPtr
)
{
    if (mesh.debug)
    {
        InfoInFunction << "Checking points" << endl;
    }

    label nFaceErrors = 0;
    label nCellErrors = 0;

    const labelListList& pf = mesh.pointFaces();

    forAll(pf, pointi)
    {
        if (pf[pointi].empty())
        {
            if (setPtr)
            {
                setPtr->insert(pointi);
            }

            nFaceErrors++;
        }
    }


    forAll(pf, pointi)
    {
        const labelList& pc = mesh.pointCells(pointi);

        if (pc.empty())
        {
            if (setPtr)
            {
                setPtr->insert(pointi);
            }

            nCellErrors++;
        }
    }

    reduce(nFaceErrors, sumOp<label>());
    reduce(nCellErrors, sumOp<label>());

    if (nFaceErrors > 0 || nCellErrors > 0)
    {
        if (report)
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
        if (report)
        {
            Info<< "    Point usage OK." << endl;
        }

        return false;
    }
}


bool Foam::meshCheck::checkDuplicateFaces
(
    const primitiveMesh& mesh,
    const label facei,
    const Map<label>& nCommonPoints,
    label& nBaffleFaces,
    labelHashSet* setPtr
)
{
    bool error = false;

    forAllConstIter(Map<label>, nCommonPoints, iter)
    {
        const label nbFacei = iter.key();
        const label nCommon = iter();

        const face& curFace = mesh.faces()[facei];
        const face& nbFace = mesh.faces()[nbFacei];

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
                setPtr->insert(facei);
                setPtr->insert(nbFacei);
            }
        }
    }

    return error;
}


bool Foam::meshCheck::checkCommonOrder
(
    const primitiveMesh& mesh,
    const label facei,
    const Map<label>& nCommonPoints,
    labelHashSet* setPtr
)
{
    bool error = false;

    forAllConstIter(Map<label>, nCommonPoints, iter)
    {
        const label nbFacei = iter.key();
        const label nCommon = iter();

        const face& curFace = mesh.faces()[facei];
        const face& nbFace = mesh.faces()[nbFacei];

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
                    const label fpPlus1 = curFace.fcIndex(fp);
                    const label fpMin1  = curFace.rcIndex(fp);

                    // Vertices before and after on nbFace
                    const label nbPlus1 = nbFace.fcIndex(nb);
                    const label nbMin1  = nbFace.rcIndex(nb);

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
                                setPtr->insert(facei);
                                setPtr->insert(nbFacei);
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


bool Foam::meshCheck::checkFaceFaces
(
    const primitiveMesh& mesh,
    const bool report,
    labelHashSet* setPtr
)
{
    if (mesh.debug)
    {
        InfoInFunction << "Checking face-face connectivity" << endl;
    }

    const labelListList& pf = mesh.pointFaces();

    label nBaffleFaces = 0;
    label nErrorDuplicate = 0;
    label nErrorOrder = 0;
    Map<label> nCommonPoints(100);

    for (label facei = 0; facei < mesh.nFaces(); facei++)
    {
        const face& curFace = mesh.faces()[facei];

        // Calculate number of common points between current facei and
        // neighbouring face. Store on map.
        nCommonPoints.clear();

        forAll(curFace, fp)
        {
            const label pointi = curFace[fp];
            const labelList& nbs = pf[pointi];

            forAll(nbs, nbI)
            {
                const label nbFacei = nbs[nbI];

                if (facei < nbFacei)
                {
                    // Only check once for each combination of two faces.

                    Map<label>::iterator fnd = nCommonPoints.find(nbFacei);

                    if (fnd == nCommonPoints.end())
                    {
                        // First common vertex found.
                        nCommonPoints.insert(nbFacei, 1);
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
        if
        (
            checkDuplicateFaces
            (
                mesh,
                facei,
                nCommonPoints,
                nBaffleFaces,
                setPtr
            )
        )
        {
            nErrorDuplicate++;
        }

        // Check common vertices are consecutive on both faces
        if (checkCommonOrder(mesh, facei, nCommonPoints, setPtr))
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

        return false;   // return true;
    }
    else
    {
        if (report)
        {
            Info<< "    Face-face connectivity OK." << endl;
        }

        return false;
    }
}


// ************************************************************************* //
