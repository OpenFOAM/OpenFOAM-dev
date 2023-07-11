/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2020-2023 OpenFOAM Foundation
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

#include "MPLICcell.H"
#include "tetCell.H"
#include "cubicEqn.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::label Foam::MPLICcell::calcMatchAlphaCutCell
(
    const MPLICcellStorage& cellInfo,
    const bool tetDecom
)
{
    // Clear all the temporary lists/fields
    clear();

    const scalar cellAlpha = cellInfo.cellAlpha();

    const UIndirectList<scalar>& cellMagSfs = cellInfo.magSf();

    // Initialise fluxes from velocity point values
    if (!unweighted_)
    {
        phiU
        (
            cellInfo.points(),
            cellInfo.faces(),
            cellInfo.cellFaces(),
            cellInfo.pointsU()
        );
    }

    // Difference between point alpha values isn't large enough
    if (mag(cellInfo.cellAlphaMax() - cellInfo.cellAlphaMin()) < vSmall)
    {
        return -1;
    }

    // Finding between which two points cell cut lies
    findPointAlphaBounds(cellInfo, tetDecom);

    // Direction of the cut isn't clear => use default scheme
    if (mag(pCubicAlphas_.a() - pCubicAlphas_.d()) < vSmall)
    {
        return -1;
    }

    // Calculates 2 intermediate volume fractions necessary for cubic polyfit
    calcPointAlphaInterior(cellInfo, tetDecom);

    // Calculate coefficients of cubic polynomial (computed on normalised input)
    const FixedList<scalar, 4> coeffs(solveVanderMatrix());

    // Direct solution of cubic equation and roots selection
    findRoots(cellInfo, coeffs, tetDecom);

    // Use default scheme if the calculated alpha is more than 10% off from
    // original cell alpha.
    if (mag(cutAlpha_) > rootSmall && (1 - mag(cellAlpha/cutAlpha_)) > 0.1)
    {
        return 0;
    }
    else
    {
        unweighted_ ? calcAlphaf(cellMagSfs) : calcAlphaUf();
        return 1;
    }
}


void Foam::MPLICcell::findPointAlphaBounds
(
    const MPLICcellStorage& cellInfo,
    const bool tetDecom
)
{
    const scalar cellAlpha = cellInfo.cellAlpha();

    cellPointsAlpha_ =
        UIndirectList<scalar>(cellInfo.pointsAlpha(), cellInfo.cellPoints());

    if (tetDecom)
    {
        cellPointsAlpha_.append(cellAlpha);
    }

    sort(cellPointsAlpha_);

    // Avoid useless cuts and keep only unique point values

    pointsAlpha_.clear();
    pointsAlpha_.append(cellPointsAlpha_[0]);
    for (label i=1; i<cellPointsAlpha_.size(); i++)
    {
        if (mag(cellPointsAlpha_[i-1] - cellPointsAlpha_[i]) > vSmall)
        {
            pointsAlpha_.append(cellPointsAlpha_[i]);
        }
    }

    const label nAlphas = pointsAlpha_.size();

    // If there is only one value (all the point are the same)
    // the cut is not determined
    if (nAlphas < 2)
    {
        pCubicAlphas_.a() = -1;
        pCubicAlphas_.d() = -1;
        cCubicAlphas_.a() = -1;
        cCubicAlphas_.d() = -1;
    }

    // If in the cell are at least two points with different values we can
    // attempt to cut between them
    else if (nAlphas == 2)
    {
        pCubicAlphas_.a() = pointsAlpha_.first();
        pCubicAlphas_.d() = pointsAlpha_.last();
        cCubicAlphas_.a() = 1;
        cCubicAlphas_.d() = 0;
    }

    // We know extreme of the volume fraction so for the points that have those
    // it will be only 0 and 1, therefore we make only necessary cuts starting
    // from the mid index point value
    else
    {
        // Pick the mid point value
        label index = label(round((nAlphas)/2.0)) - 1;

        // Calculate initial cell and point value
        scalar target = pointsAlpha_[index];
        scalar cutAlpha = calcAlpha(cellInfo, target, tetDecom);
        scalar targetOld = target;
        scalar cutAlphaOld = cutAlpha;

        for (label i = 1; i < nAlphas-1; ++i)
        {
            (cutAlpha >= cellAlpha) ? index++ : index--;

            // Special case
            // Maximum point value and minimum volume fraction = 0
            if (index == nAlphas - 1)
            {
                cCubicAlphas_.a() = cutAlpha;
                cCubicAlphas_.d() = 0;
                pCubicAlphas_.a() = target;
                pCubicAlphas_.d() = pointsAlpha_[index];
                break;
            }

            // Special case
            // Minimum point value and maximum volume fraction = 1
            else if (index == 0)
            {
                cCubicAlphas_.a() = 1;
                cCubicAlphas_.d() = cutAlpha;
                pCubicAlphas_.a() = pointsAlpha_[index];
                pCubicAlphas_.d() = target;
                break;
            }

            // Calculate new values
            target = pointsAlpha_[index];
            cutAlpha = calcAlpha(cellInfo, target, tetDecom);

            if (cutAlphaOld > cellAlpha && cutAlpha < cellAlpha)
            {
                cCubicAlphas_.a() = cutAlphaOld;
                cCubicAlphas_.d() = cutAlpha;
                pCubicAlphas_.a() = targetOld;
                pCubicAlphas_.d() = target;
                break;
            }

            else if (cutAlphaOld < cellAlpha && cutAlpha > cellAlpha)
            {
                cCubicAlphas_.a() = cutAlpha;
                cCubicAlphas_.d() = cutAlphaOld;
                pCubicAlphas_.a() = target;
                pCubicAlphas_.d() = targetOld;
                break;
            }

            // Store previous iteration values
            targetOld = target;
            cutAlphaOld = cutAlpha;
        }
    }
}


void Foam::MPLICcell::calcPointAlphaInterior
(
    const MPLICcellStorage& cellInfo,
    const bool tetDecom
)
{
    for (label i=1; i<=2; i++)
    {
        pCubicAlphas_[i] =
            pCubicAlphas_.a() + (pCubicAlphas_.d() - pCubicAlphas_.a())*(i/3.0);

        cCubicAlphas_[i] = calcAlpha(cellInfo, pCubicAlphas_[i], tetDecom);
    }
}


Foam::FixedList<Foam::scalar, 4> Foam::MPLICcell::solveVanderMatrix() const
{
    // The cubic polynomial of the volume is fit to a normalised coordinate,
    // which is defined as follows (see the use of sFactor in findRoots below):
    //
    //      x = (alpha - pCubicAlpha[0])/(pCubicAlpha[3] - pCubicAlpha[0])
    //
    // pCubicAlpha[0] and pCubicAlpha[3] are the bounds of the fit, and are
    // actual iso-values on points of the cell. pCubicAlpha_1 and pCubicAlpha_2
    // are interior values needed to complete the cubic fit. The 4 values
    // are equidistant (see calcPointAlphaInterior above) so the corresponding
    // values of the normalised coordinate x are always:
    //
    //      x = [0   1/3  2/3 1]
    //
    // This means the Vandermonde matrix that is solved for the cubic
    // polynomial's coefficients is always the same:
    //
    //      V = [0    0   0   1]
    //          [1/27 1/9 1/3 1]
    //          [8/27 4/9 2/3 1]
    //          [1    1   1   1]
    //
    // This means its inverse can be precomputed. This pre-computation is hard
    // coded below.

    const vector4& b = cCubicAlphas_;

    return FixedList<scalar, 4>
    {
        scalar(- 4.5*b[0] + 13.5*b[1] - 13.5*b[2] + 4.5*b[3]),
        scalar(9.0*b[0] - 22.5*b[1] + 18.0*b[2] - 4.5*b[3]),
        scalar(- 5.5*b[0] +  9.0*b[1] -  4.5*b[2] + 1.0*b[3]),
        scalar(1.0*b[0])
    };
}


void Foam::MPLICcell::findRoots
(
    const MPLICcellStorage& cellInfo,
    const FixedList<scalar, 4>& coeff,
    const bool tetDecom
)
{
    const scalar cellAlpha = cellInfo.cellAlpha();

    // Solve cubic polynomial exactly
    const Roots<3> roots =
        cubicEqn(coeff[0], coeff[1], coeff[2], coeff[3] - cellAlpha).roots();

    // Find which root corresponds to the desired value
    scalar rootOld = SMALL;
    scalar target = 0;

    label nRoots = 0;
    Roots<3> selectedRoots;

    const scalar pMax = cmptMax(pCubicAlphas_);
    const scalar pMin = cmptMin(pCubicAlphas_);
    const scalar sFactor = pCubicAlphas_.d() - pCubicAlphas_.a();

    forAll(roots, rooti)
    {
        // Scale the roots back into original scale
        const scalar root = roots[rooti]*sFactor + pCubicAlphas_.a();

        // Pick up correct root
        if (root < pMax && root > pMin && rootOld != root)
        {
            target = (target == 0) ? root : target;
            selectedRoots[nRoots++] = root;
        }

        rootOld = root;
    }

    // Recompute alpha last for analytical approach
    cutAlpha_ = calcAlpha(cellInfo, target, tetDecom);

    // In case the selection of the root failed compute all three volume
    // fractions and choose the one with minimum error
    scalar error = mag(cutAlpha_ - cellAlpha);

    // If error > 1e-3 check for better root
    if (nRoots > 0 && error > 1e-3)
    {
        scalar minError = error;
        label minIndex = 0;

        for (label rooti=1; rooti<nRoots; rooti++)
        {
            const scalar targeti =
                calcAlpha(cellInfo, selectedRoots[rooti], tetDecom);

            error = mag(cellAlpha - targeti);

            if (error < minError)
            {
                minError = error;
                minIndex = rooti;
            }
        }

        cutAlpha_ = calcAlpha(cellInfo, selectedRoots[minIndex], tetDecom);
    }
}


Foam::scalar Foam::MPLICcell::calcAlpha
(
    const MPLICcellStorage& cellInfo,
    const scalar target,
    const bool tetDecom
)
{
    if (!tetDecom)
    {
        return calcCutCellVolumeAlpha(cellInfo, target);
    }
    else
    {
        return calcTetCutCellVolumeAlpha(cellInfo, target);
    }
}


void Foam::MPLICcell::calcSubCellVolume()
{
    vector cEst = subFaceCentres_[0];
    for(label i = 1; i < subFaceCentres_.size(); i++)
    {
        cEst += subFaceCentres_[i];
    }
    cEst /= subFaceCentres_.size();

    subCellVolume_ = 0;
    forAll(subFaceAreas_, i)
    {
        subCellVolume_ += subFaceAreas_[i] & (subFaceCentres_[i] - cEst);
    }
    subCellVolume_ /= 3.0;
}


Foam::scalar Foam::MPLICcell::calcCutCellVolumeAlpha
(
    const MPLICcellStorage& cellInfo,
    const scalar target
)
{
    const scalar V = cellInfo.V();

    // Case when the cell has to be cut
    if (cellInfo.cellAlphaMax() > target && cellInfo.cellAlphaMin() < target)
    {
        // Cut cell single cut if multicut detected use multicut
        const bool status = singleCutCell(cellInfo, target);
        if (!status && multiCut_)
        {
            multiCutCell(cellInfo, target);
        }

        // Compute normal
        cutNormal_ = normalised(cutSf_);

        // Calculate volume
        if (subFaceCentres_.size() != 0)
        {
            calcSubCellVolume();
        }

        // Snap negative volume cell to zero
        if (subCellVolume_ <= 0)
        {
            resetFaceFields(cellInfo.size());
            subCellVolume_ = 0;
            return 0;
        }

        return min(subCellVolume_, V)/V;
    }
    else if (target <= cellInfo.cellAlphaMin())
    {
        if (unweighted_)
        {
            subFaceMagSf_ = cellInfo.magSf();
        }
        else
        {
            alphaPhiU_ = phiU_;
        }
        subCellVolume_ = V;

        return 1;
    }
    else
    {
        resetFaceFields(cellInfo.size());
        subCellVolume_ = 0;

        return 0;
    }
}


Foam::scalar Foam::MPLICcell::calcTetCutCellVolumeAlpha
(
    const MPLICcellStorage& cellInfo,
    const scalar target
)
{
    clear();
    resetFaceFields(cellInfo.size());

    // Append cell centre value
    pointsAlpha_ = cellInfo.pointsAlpha();
    pointsAlpha_.append(target);

    // Overall volume
    scalar cellVolume = 0;

    if (min(pointsAlpha_) < target && max(pointsAlpha_) > target)
    {
        // Cell centre is the first point of the tet for all tets in the cell
        const vector& a = cellInfo.C();

        // Cell centre value is the first value of the tet
        // for all tets in the cell
        tetPointsAlpha_[0] = cellInfo.cellAlpha();

        if (!unweighted_)
        {
            tetPointsU_[0] = cellInfo.cellU();
        }

        // Looping through all the faces
        forAll(cellInfo, facei)
        {
            // Create copy of the face indexing in order to flip if necessary
            face f = cellInfo.faces()[cellInfo.cellFaces()[facei]];

            // Work directly with all the faces pointing out of the cell
            if (!cellInfo.isOwner()[facei])
            {
                f.flip();
            }

            const label& bL = f[0];
            const point& b = cellInfo.points()[bL];
            tetPointsAlpha_[1] = cellInfo.pointsAlpha()[bL];
            if (!unweighted_)
            {
                tetPointsU_[1] = cellInfo.pointsU()[bL];
            }

            // Decomposing faces
            for (label i = 1; i < f.size()-1; ++i)
            {
                // Labels for point c and d
                const label cL = f[i];
                const label dL = f[i + 1];

                // c, d points of tetrahedron
                const point& c = cellInfo.points()[cL];
                const point& d = cellInfo.points()[dL];

                // Tet point values
                tetPointsAlpha_[2] = cellInfo.pointsAlpha()[cL];
                tetPointsAlpha_[3] = cellInfo.pointsAlpha()[dL];

                if (!unweighted_)
                {
                    tetPointsU_[2] = cellInfo.pointsU()[cL];
                    tetPointsU_[3] = cellInfo.pointsU()[dL];
                }

                // Tet maximum and minimum point values
                const scalar tetMax = max(tetPointsAlpha_);
                const scalar tetMin = min(tetPointsAlpha_);

                // Contains all the geometric information
                tetPointRef cellTet(a, b, c, d);

                // Integrate overall volume for consistency
                cellVolume += cellTet.mag();

                // Tet cuts
                if (tetMin < target && tetMax > target)
                {
                    // Tet point
                    tetPoints_ = {a, b, c, d};
                    tetSf_[0] = cellTet.Sa();
                    tetSf_[1] = cellTet.Sb();
                    tetSf_[2] = cellTet.Sc();
                    tetSf_[3] = cellTet.Sd();
                    tetCf_[0] = triPointRef(b, c, d).centre();
                    tetCf_[1] = triPointRef(a, d, c).centre();
                    tetCf_[2] = triPointRef(a, b, d).centre();
                    tetCf_[3] = triPointRef(a, c, b).centre();

                    // Geometric cut of tetrahedra cell
                    const bool ow = cellInfo.isOwner()[facei];

                    // Tetrahedron cut
                    cutTetCell(target, facei,ow);

                    if (subFaceCentres_.size() > 0)
                    {
                        calcSubCellVolume();
                    }
                }

                // Fully submerged tet
                else if (tetMin >= target)
                {
                    subCellVolume_ += cellTet.mag();
                    if (unweighted_)
                    {
                        subFaceMagSf_[facei] += mag(cellTet.Sa());
                    }
                    else
                    {
                        const scalar phiU =
                            (
                                (1.0/3.0)*
                                (
                                    tetPointsU_[1]
                                  + tetPointsU_[2]
                                  + tetPointsU_[3]
                                )
                            ) & cellTet.Sa();

                        if (cellInfo.isOwner()[facei])
                        {
                            alphaPhiU_[facei] += phiU;
                        }
                        else
                        {
                            alphaPhiU_[facei] -= phiU;
                        }
                    }
                }
            }
        }

        // Compute normal
        cutNormal_ = normalised(cutSf_);

        // Snap negative volume cell to zero
        if (subCellVolume_ <= 0)
        {
            resetFaceFields(cellInfo.size());
            subCellVolume_ = 0;
            return 0;
        }
        return min(subCellVolume_, cellVolume)/cellVolume;
    }
    else if (target <= min(pointsAlpha_))
    {
        if (unweighted_)
        {
            subFaceMagSf_ = cellInfo.magSf();
        }
        else
        {
            alphaPhiU_ = phiU_;
        }
        subCellVolume_ = cellVolume;
        return 1;
    }
    else
    {
        resetFaceFields(cellInfo.size());
        subCellVolume_ = 0;
        return 0;
    }
}


bool Foam::MPLICcell::singleCutCell
(
    const MPLICcellStorage& cellInfo,
    const scalar target
)
{
    clear();
    resetFaceFields(cellInfo.size());

    // Cut type
    label cutType;

    // Any face has more then one cut?
    bool moreCutsPerFace = 0;

    // Single cell cut
    forAll(cellInfo, facei)
    {
        // Collect fully submerged faces
        if (cellInfo.facesAlphaMin()[facei] >= target)
        {
            appendSfCf
            (
                cellInfo.Sf()[facei],
                cellInfo.Cf()[facei],
                cellInfo.magSf()[facei],
                cellInfo.isOwner()[facei]
            );

            if (unweighted_)
            {
                subFaceMagSf_[facei] = cellInfo.magSf()[facei];
            }
            else
            {
                alphaPhiU_[facei] = phiU_[facei];
            }
            continue;
        }
        else if (cellInfo.facesAlphaMax()[facei] < target)
        {
            continue;
        }

        // Cut the face return label of next face and edge
        cutType = faceCutter_.cutFace
        (
            cellInfo.faces()[cellInfo.cellFaces()[facei]],
            cellInfo.points(),
            cellInfo.pointsAlpha(),
            cellInfo.pointsU(),
            target,
            cellInfo.isOwner()[facei]
        );

        // Potentially multiple cuts through the cell
        if (cutType == -1)
        {
            moreCutsPerFace = 1;
        }

        else if (cutType == 1)
        {
            // Append to the cut list of points
            cutPoints_.append(faceCutter_.cutPoints());

            // Append area vectors and face centers
            if (faceCutter_.subPoints().size() > 2)
            {
                const vector Sf = faceCutter_.Sf();
                const vector Cf = faceCutter_.Cf(Sf);
                const scalar magSf = mag(Sf);
                appendSfCf(Sf, Cf, magSf);

                if (unweighted_)
                {
                    subFaceMagSf_[facei] += magSf;
                }
                else
                {
                    alphaPhiU_[facei] += faceCutter_.alphaPhiU();
                }
            }
        }
    }

    // Assume it is multicut if triangle have opposite sign in any direction
    bool cutOrientationDiffers = 0;
    if (cutPoints_.size() > 2)
    {
        cutOrientationDiffers = cutStatusCalcSf();
        const vector Cf = calcCutCf(cutSf_);
        appendSfCf(cutSf_, Cf, mag(cutSf_));
    }

    // Potentially multiple cuts through cell
    if (cutOrientationDiffers || moreCutsPerFace)
    {
        return 0;
    }

    // Only one cut through the cell
    else
    {
        return 1;
    }
}


bool Foam::MPLICcell::multiCutCell
(
    const MPLICcellStorage& cellInfo,
    const scalar target
)
{
    clear();
    resetFaceFields(cellInfo.size());

    // Prepare local addressing
    if (!addressingCalculated_)
    {
        calcAddressing(cellInfo);
    }

    // Keep track of cut edges
    boolList isEdgeCutOld(cellInfo.cellEdges().size(), false);
    boolList isEdgeCut(cellInfo.cellEdges().size(), false);

    // Keep track of the fully submerged subfaces
    boolList submerged(cellInfo.size(), false);

    // Initialise the list of necessary labels
    label facei, nextFace, faceEdgei, status;

    // Loop through all the cuts
    // Number of cuts limited to number of faces
    forAll(cellInfo, cutI)
    {
        faceEdgei = -1;
        facei = 0;
        nextFace = 0;
        status = 0;

        // One cell cut
        label j = 0;
        while (j < cellInfo.size())
        {
            facei = (status == 0) ? j : nextFace;

            // Collect fully submerged faces
            if (cellInfo.facesAlphaMin()[facei] >= target && !submerged[facei])
            {
                submerged[facei] = true;
                appendSfCf
                (
                    cellInfo.Sf()[facei],
                    cellInfo.Cf()[facei],
                    cellInfo.magSf()[facei],
                    cellInfo.isOwner()[facei]
                );

                // Precompute face fields
                if (unweighted_)
                {
                    subFaceMagSf_[facei] = cellInfo.magSf()[facei];
                }
                else
                {
                    alphaPhiU_[facei] = phiU_[facei];
                }
            }

            // Cut the face
            status = faceCutter_.cutFace
            (
                cellInfo.faces()[cellInfo.cellFaces()[facei]],
                localFaceEdges_[facei],
                cellInfo.points(),
                isEdgeCutOld,
                isEdgeCut,
                faceEdgei,
                cellInfo.pointsAlpha(),
                cellInfo.pointsU(),
                facei,
                target,
                cellInfo.isOwner()[facei]
            );

            // Get the next face and edge
            if (status)
            {
                const label edgei = localFaceEdges_[facei][faceEdgei];
                const labelList& edgeFaces = localEdgeFaces_[edgei];
                nextFace = edgeFaces[edgeFaces[0] == facei];
                faceEdgei = findIndex(localFaceEdges_[nextFace], edgei);
            }

            // Append to the cut list of points
            cutPoints_.append(faceCutter_.cutPoints());
            cutEdges_.append(faceCutter_.cutEdges());

            // Append area vectors and face centers
            if (faceCutter_.subPoints().size() > 2 && !submerged[facei])
            {
                const vector Sf = faceCutter_.Sf();
                const vector Cf = faceCutter_.Cf(Sf);
                const scalar magSf = mag(Sf);

                // The sub-faces are always pointing outwards
                appendSfCf(Sf, Cf, magSf);

                if (unweighted_)
                {
                    subFaceMagSf_[facei] += magSf;
                }
                else
                {
                    alphaPhiU_[facei] += faceCutter_.alphaPhiU();
                }
            }

            // End on reaching first edge
            if (cutEdges_.size() > 0 && cutEdges_.first() == cutEdges_.last())
            {
                break;
            }
            if (status == 0)
            {
                ++j;
            }
        }

        isEdgeCutOld = isEdgeCut;

        if (cutPoints_.size() == 0)
        {
            break;
        }
        else
        {
            // Append information from cut face
            const vector Sf = calcCutSf();
            const vector Cf = calcCutCf(Sf);
            appendSfCf(Sf, Cf, mag(Sf));
            cutSf_ += Sf;
        }

        // Clear fields to prepare for next cut
        cutPoints_.clear();
        cutEdges_.clear();
    }
    return 1;
}


bool Foam::MPLICcell::cutTetCell
(
    const scalar target,
    const label faceOrig,
    const bool ow
)
{
    // Clear geometry data for tet cut
    clearOneCut();

    // Single cell cut
    forAll(tetFaces_, facei)
    {
        const face& f = tetFaces_[facei];

        // Collect fully submerged faces
        if
        (
            min
            (
                min
                (
                    tetPointsAlpha_[f[0]],
                    tetPointsAlpha_[f[1]]
                ),
                tetPointsAlpha_[f[2]]
            ) >= target
        )
        {
            const vector& Sf = tetSf_[facei];
            const vector& Cf = tetCf_[facei];
            appendSfCf(Sf, Cf, mag(Sf));

            if (unweighted_ && facei == 0)
            {
                subFaceMagSf_[faceOrig] += mag(Sf);
            }
            else if (!unweighted_ && facei == 0)
            {
                const face& f0 = tetFaces_[0];
                const scalar phiU =
                    (
                        (1.0/3.0)*
                        (
                            tetPointsU_[f0[0]]
                          + tetPointsU_[f0[1]]
                          + tetPointsU_[f0[2]]
                        )
                    ) & Sf;
                if (ow)
                {
                    alphaPhiU_[faceOrig] += phiU;
                }
                else
                {
                    alphaPhiU_[faceOrig] -= phiU;
                }
            }
            continue;
        }
        else if
        (
            max
            (
                max
                (
                    tetPointsAlpha_[f[0]],
                    tetPointsAlpha_[f[1]]
                ),
                tetPointsAlpha_[f[2]]
            ) < target
        )
        {
            continue;
        }

        // Cut the face return label of next face and edge
        faceCutter_.cutFace
        (
            tetFaces_[facei],
            tetPoints_,
            tetPointsAlpha_,
            tetPointsU_,
            target,
            true
        );

        // Append to the cut list of points
        cutPoints_.append(faceCutter_.cutPoints());

        // Append area vectors and face centers
        if (faceCutter_.subPoints().size() > 2)
        {
            const vector Sf = faceCutter_.Sf();
            const vector Cf = faceCutter_.Cf(Sf);
            const scalar magSf = mag(Sf);
            appendSfCf(Sf, Cf, magSf);

            // For unweighted alphaf
            if (unweighted_ && facei == 0)
            {
                subFaceMagSf_[faceOrig] += magSf;
            }

            // For phiU weighted alphaf
            else if (!unweighted_ && facei == 0)
            {
                if (ow)
                {
                    alphaPhiU_[faceOrig] += faceCutter_.alphaPhiU();
                }
                else
                {
                    alphaPhiU_[faceOrig] -= faceCutter_.alphaPhiU();
                }
            }
        }
    }

    // Append information from cut face
    if (cutPoints_.size() > 2)
    {
        const vector Sf = calcCutSf();
        const vector Cf = calcCutCf(Sf);
        appendSfCf(Sf, Cf, mag(Sf));
        cutSf_  += Sf;
    }

    return 1;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::MPLICcell::MPLICcell(const bool unweighted, const bool multiCut)
:
    unweighted_(unweighted),
    multiCut_(multiCut),
    faceCutter_(unweighted),
    cutPoints_(10),
    cutEdges_(10),
    subFaceAreas_(10),
    subFaceCentres_(10),
    tetPointsAlpha_(4),
    tetPointsU_(4),
    tetFaces_
    {
        triFace(1, 2, 3),
        triFace(0, 3, 2),
        triFace(0, 1, 3),
        triFace(0, 2, 1)
    },
    pointsAlpha_(8)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::MPLICcell::matchAlpha
(
    const MPLICcellStorage& cellInfo
)
{
    // Addressing for multicut needs to be recomputed for each cell
    addressingCalculated_ = false;

    // Try normal cell cut matching first
    label status = calcMatchAlphaCutCell(cellInfo);

    // If volume fraction error is bigger than 10% try tet decomposition cut
    if (status == 0 && multiCut_)
    {
        status = calcMatchAlphaCutCell(cellInfo, true);
    }

    if (status == 0 || status == -1)
    {
        return 0;
    }
    else
    {
        return 1;
    }
}


// ************************************************************************* //
