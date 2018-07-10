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

\*---------------------------------------------------------------------------*/

#include "polyMeshGeometry.H"
#include "polyMeshTetDecomposition.H"
#include "pyramidPointFaceRef.H"
#include "tetPointRef.H"
#include "syncTools.H"
#include "unitConversion.H"
#include "primitiveMeshTools.H"

namespace Foam
{
    defineTypeNameAndDebug(polyMeshGeometry, 0);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::polyMeshGeometry::updateFaceCentresAndAreas
(
    const pointField& p,
    const labelList& changedFaces
)
{
    const faceList& fs = mesh_.faces();

    forAll(changedFaces, i)
    {
        label facei = changedFaces[i];

        const labelList& f = fs[facei];
        label nPoints = f.size();

        // If the face is a triangle, do a direct calculation for efficiency
        // and to avoid round-off error-related problems
        if (nPoints == 3)
        {
            faceCentres_[facei] = (1.0/3.0)*(p[f[0]] + p[f[1]] + p[f[2]]);
            faceAreas_[facei] = 0.5*((p[f[1]] - p[f[0]])^(p[f[2]] - p[f[0]]));
        }
        else
        {
            vector sumN = Zero;
            scalar sumA = 0.0;
            vector sumAc = Zero;

            point fCentre = p[f[0]];
            for (label pi = 1; pi < nPoints; pi++)
            {
                fCentre += p[f[pi]];
            }

            fCentre /= nPoints;

            for (label pi = 0; pi < nPoints; pi++)
            {
                const point& nextPoint = p[f[(pi + 1) % nPoints]];

                vector c = p[f[pi]] + nextPoint + fCentre;
                vector n = (nextPoint - p[f[pi]])^(fCentre - p[f[pi]]);
                scalar a = mag(n);

                sumN += n;
                sumA += a;
                sumAc += a*c;
            }

            faceCentres_[facei] = (1.0/3.0)*sumAc/(sumA + vSmall);
            faceAreas_[facei] = 0.5*sumN;
        }
    }
}


void Foam::polyMeshGeometry::updateCellCentresAndVols
(
    const labelList& changedCells,
    const labelList& changedFaces
)
{
    // Clear the fields for accumulation
    UIndirectList<vector>(cellCentres_, changedCells) = Zero;
    UIndirectList<scalar>(cellVolumes_, changedCells) = 0.0;

    const labelList& own = mesh_.faceOwner();
    const labelList& nei = mesh_.faceNeighbour();

    // first estimate the approximate cell centre as the average of face centres

    vectorField cEst(mesh_.nCells());
    UIndirectList<vector>(cEst, changedCells) = Zero;
    scalarField nCellFaces(mesh_.nCells());
    UIndirectList<scalar>(nCellFaces, changedCells) = 0.0;

    forAll(changedFaces, i)
    {
        label facei = changedFaces[i];
        cEst[own[facei]] += faceCentres_[facei];
        nCellFaces[own[facei]] += 1;

        if (mesh_.isInternalFace(facei))
        {
            cEst[nei[facei]] += faceCentres_[facei];
            nCellFaces[nei[facei]] += 1;
        }
    }

    forAll(changedCells, i)
    {
        label celli = changedCells[i];
        cEst[celli] /= nCellFaces[celli];
    }

    forAll(changedFaces, i)
    {
        label facei = changedFaces[i];

        // Calculate 3*face-pyramid volume
        scalar pyr3Vol = max
        (
            faceAreas_[facei] & (faceCentres_[facei] - cEst[own[facei]]),
            vSmall
        );

        // Calculate face-pyramid centre
        vector pc = (3.0/4.0)*faceCentres_[facei] + (1.0/4.0)*cEst[own[facei]];

        // Accumulate volume-weighted face-pyramid centre
        cellCentres_[own[facei]] += pyr3Vol*pc;

        // Accumulate face-pyramid volume
        cellVolumes_[own[facei]] += pyr3Vol;

        if (mesh_.isInternalFace(facei))
        {
            // Calculate 3*face-pyramid volume
            scalar pyr3Vol = max
            (
                faceAreas_[facei] & (cEst[nei[facei]] - faceCentres_[facei]),
                vSmall
            );

            // Calculate face-pyramid centre
            vector pc =
                (3.0/4.0)*faceCentres_[facei]
              + (1.0/4.0)*cEst[nei[facei]];

            // Accumulate volume-weighted face-pyramid centre
            cellCentres_[nei[facei]] += pyr3Vol*pc;

            // Accumulate face-pyramid volume
            cellVolumes_[nei[facei]] += pyr3Vol;
        }
    }

    forAll(changedCells, i)
    {
        label celli = changedCells[i];

        cellCentres_[celli] /= cellVolumes_[celli] + vSmall;
        cellVolumes_[celli] *= (1.0/3.0);
    }
}


Foam::labelList Foam::polyMeshGeometry::affectedCells
(
    const polyMesh& mesh,
    const labelList& changedFaces
)
{
    const labelList& own = mesh.faceOwner();
    const labelList& nei = mesh.faceNeighbour();

    labelHashSet affectedCells(2*changedFaces.size());

    forAll(changedFaces, i)
    {
        label facei = changedFaces[i];

        affectedCells.insert(own[facei]);

        if (mesh.isInternalFace(facei))
        {
            affectedCells.insert(nei[facei]);
        }
    }
    return affectedCells.toc();
}


Foam::scalar Foam::polyMeshGeometry::checkNonOrtho
(
    const polyMesh& mesh,
    const bool report,
    const scalar severeNonorthogonalityThreshold,
    const label facei,
    const vector& s,    // face area vector
    const vector& d,    // cc-cc vector

    label& severeNonOrth,
    label& errorNonOrth,
    labelHashSet* setPtr
)
{
    scalar dDotS = (d & s)/(mag(d)*mag(s) + vSmall);

    if (dDotS < severeNonorthogonalityThreshold)
    {
        label nei = -1;

        if (mesh.isInternalFace(facei))
        {
            nei = mesh.faceNeighbour()[facei];
        }

        if (dDotS > small)
        {
            if (report)
            {
                // Severe non-orthogonality but mesh still OK
                Pout<< "Severe non-orthogonality for face " << facei
                    << " between cells " << mesh.faceOwner()[facei]
                    << " and " << nei
                    << ": Angle = "
                    << radToDeg(::acos(dDotS))
                    << " deg." << endl;
            }

            severeNonOrth++;
        }
        else
        {
            // Non-orthogonality greater than 90 deg
            if (report)
            {
                WarningInFunction
                    << "Severe non-orthogonality detected for face "
                    << facei
                    << " between cells " << mesh.faceOwner()[facei]
                    << " and " << nei
                    << ": Angle = "
                    << radToDeg(::acos(dDotS))
                    << " deg." << endl;
            }

            errorNonOrth++;
        }

        if (setPtr)
        {
            setPtr->insert(facei);
        }
    }
    return dDotS;
}


// Create the neighbour pyramid - it will have positive volume
bool Foam::polyMeshGeometry::checkFaceTet
(
    const polyMesh& mesh,
    const bool report,
    const scalar minTetQuality,
    const pointField& p,
    const label facei,
    const point& fc,    // face centre
    const point& cc,    // cell centre

    labelHashSet* setPtr
)
{
    const face& f = mesh.faces()[facei];

    forAll(f, fp)
    {
        scalar tetQual = tetPointRef
        (
            p[f[fp]],
            p[f.nextLabel(fp)],
            fc,
            cc
        ).quality();

        if (tetQual < minTetQuality)
        {
            if (report)
            {
                Pout<< "bool polyMeshGeometry::checkFaceTets("
                    << "const bool, const scalar, const pointField&"
                    << ", const pointField&"
                    << ", const labelList&, labelHashSet*) : "
                    << "face " << facei
                    << " has a triangle that points the wrong way."
                     << endl
                    << "Tet quality: " << tetQual
                    << " Face " << facei
                    << endl;
            }

            if (setPtr)
            {
                setPtr->insert(facei);
            }
            return true;
        }
    }
    return false;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Foam::polyMeshGeometry::polyMeshGeometry(const polyMesh& mesh)
:
    mesh_(mesh)
{
    correct();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::polyMeshGeometry::correct()
{
    faceAreas_ = mesh_.faceAreas();
    faceCentres_ = mesh_.faceCentres();
    cellCentres_ = mesh_.cellCentres();
    cellVolumes_ = mesh_.cellVolumes();
}


void Foam::polyMeshGeometry::correct
(
    const pointField& p,
    const labelList& changedFaces
)
{
    // Update face quantities
    updateFaceCentresAndAreas(p, changedFaces);
    // Update cell quantities from face quantities
    updateCellCentresAndVols(affectedCells(mesh_, changedFaces), changedFaces);
}


bool Foam::polyMeshGeometry::checkFaceDotProduct
(
    const bool report,
    const scalar orthWarn,
    const polyMesh& mesh,
    const vectorField& cellCentres,
    const vectorField& faceAreas,
    const labelList& checkFaces,
    const List<labelPair>& baffles,
    labelHashSet* setPtr
)
{
    // for all internal and coupled faces check theat the d dot S product
    // is positive

    const labelList& own = mesh.faceOwner();
    const labelList& nei = mesh.faceNeighbour();
    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    // Severe nonorthogonality threshold
    const scalar severeNonorthogonalityThreshold = ::cos(degToRad(orthWarn));

    // Calculate coupled cell centre
    pointField neiCc(mesh.nFaces() - mesh.nInternalFaces());

    for (label facei = mesh.nInternalFaces(); facei < mesh.nFaces(); facei++)
    {
        neiCc[facei-mesh.nInternalFaces()] = cellCentres[own[facei]];
    }

    syncTools::swapBoundaryFacePositions(mesh, neiCc);

    scalar minDDotS = great;

    scalar sumDDotS = 0;
    label nDDotS = 0;

    label severeNonOrth = 0;

    label errorNonOrth = 0;

    forAll(checkFaces, i)
    {
        label facei = checkFaces[i];

        const point& ownCc = cellCentres[own[facei]];

        if (mesh.isInternalFace(facei))
        {
            scalar dDotS = checkNonOrtho
            (
                mesh,
                report,
                severeNonorthogonalityThreshold,
                facei,
                faceAreas[facei],
                cellCentres[nei[facei]] - ownCc,

                severeNonOrth,
                errorNonOrth,
                setPtr
            );

            if (dDotS < minDDotS)
            {
                minDDotS = dDotS;
            }

            sumDDotS += dDotS;
            nDDotS++;
        }
        else
        {
            label patchi = patches.whichPatch(facei);

            if (patches[patchi].coupled())
            {
                scalar dDotS = checkNonOrtho
                (
                    mesh,
                    report,
                    severeNonorthogonalityThreshold,
                    facei,
                    faceAreas[facei],
                    neiCc[facei-mesh.nInternalFaces()] - ownCc,

                    severeNonOrth,
                    errorNonOrth,
                    setPtr
                );

                if (dDotS < minDDotS)
                {
                    minDDotS = dDotS;
                }

                sumDDotS += dDotS;
                nDDotS++;
            }
        }
    }

    forAll(baffles, i)
    {
        label face0 = baffles[i].first();
        label face1 = baffles[i].second();

        const point& ownCc = cellCentres[own[face0]];

        scalar dDotS = checkNonOrtho
        (
            mesh,
            report,
            severeNonorthogonalityThreshold,
            face0,
            faceAreas[face0],
            cellCentres[own[face1]] - ownCc,

            severeNonOrth,
            errorNonOrth,
            setPtr
         );

        if (dDotS < minDDotS)
        {
            minDDotS = dDotS;
        }

        sumDDotS += dDotS;
        nDDotS++;
    }

    reduce(minDDotS, minOp<scalar>());
    reduce(sumDDotS, sumOp<scalar>());
    reduce(nDDotS, sumOp<label>());
    reduce(severeNonOrth, sumOp<label>());
    reduce(errorNonOrth, sumOp<label>());

    // Only report if there are some internal faces
    if (nDDotS > 0)
    {
        if (report && minDDotS < severeNonorthogonalityThreshold)
        {
            Info<< "Number of non-orthogonality errors: " << errorNonOrth
                << ". Number of severely non-orthogonal faces: "
                << severeNonOrth  << "." << endl;
        }
    }

    if (report)
    {
        if (nDDotS > 0)
        {
            Info<< "Mesh non-orthogonality Max: "
                << radToDeg(::acos(minDDotS))
                << " average: " << radToDeg(::acos(sumDDotS/nDDotS))
                << endl;
        }
    }

    if (errorNonOrth > 0)
    {
        if (report)
        {
            SeriousErrorInFunction
                << "Error in non-orthogonality detected" << endl;
        }

        return true;
    }
    else
    {
        if (report)
        {
            Info<< "Non-orthogonality check OK.\n" << endl;
        }

        return false;
    }
}


bool Foam::polyMeshGeometry::checkFacePyramids
(
    const bool report,
    const scalar minPyrVol,
    const polyMesh& mesh,
    const vectorField& cellCentres,
    const pointField& p,
    const labelList& checkFaces,
    const List<labelPair>& baffles,
    labelHashSet* setPtr
)
{
    // check whether face area vector points to the cell with higher label
    const labelList& own = mesh.faceOwner();
    const labelList& nei = mesh.faceNeighbour();

    const faceList& f = mesh.faces();

    label nErrorPyrs = 0;

    forAll(checkFaces, i)
    {
        label facei = checkFaces[i];

        // Create the owner pyramid - it will have negative volume
        scalar pyrVol = pyramidPointFaceRef
        (
            f[facei],
            cellCentres[own[facei]]
        ).mag(p);

        if (pyrVol > -minPyrVol)
        {
            if (report)
            {
                Pout<< "bool polyMeshGeometry::checkFacePyramids("
                    << "const bool, const scalar, const pointField&"
                    << ", const labelList&, labelHashSet*): "
                    << "face " << facei << " points the wrong way. " << endl
                    << "Pyramid volume: " << -pyrVol
                    << " Face " << f[facei] << " area: " << f[facei].mag(p)
                    << " Owner cell: " << own[facei] << endl
                    << "Owner cell vertex labels: "
                    << mesh.cells()[own[facei]].labels(f)
                    << endl;
            }


            if (setPtr)
            {
                setPtr->insert(facei);
            }

            nErrorPyrs++;
        }

        if (mesh.isInternalFace(facei))
        {
            // Create the neighbour pyramid - it will have positive volume
            scalar pyrVol =
                pyramidPointFaceRef(f[facei], cellCentres[nei[facei]]).mag(p);

            if (pyrVol < minPyrVol)
            {
                if (report)
                {
                    Pout<< "bool polyMeshGeometry::checkFacePyramids("
                        << "const bool, const scalar, const pointField&"
                        << ", const labelList&, labelHashSet*): "
                        << "face " << facei << " points the wrong way. " << endl
                        << "Pyramid volume: " << -pyrVol
                        << " Face " << f[facei] << " area: " << f[facei].mag(p)
                        << " Neighbour cell: " << nei[facei] << endl
                        << "Neighbour cell vertex labels: "
                        << mesh.cells()[nei[facei]].labels(f)
                        << endl;
                }

                if (setPtr)
                {
                    setPtr->insert(facei);
                }

                nErrorPyrs++;
            }
        }
    }

    forAll(baffles, i)
    {
        label face0 = baffles[i].first();
        label face1 = baffles[i].second();

        const point& ownCc = cellCentres[own[face0]];

       // Create the owner pyramid - it will have negative volume
        scalar pyrVolOwn = pyramidPointFaceRef
        (
            f[face0],
            ownCc
        ).mag(p);

        if (pyrVolOwn > -minPyrVol)
        {
            if (report)
            {
                Pout<< "bool polyMeshGeometry::checkFacePyramids("
                    << "const bool, const scalar, const pointField&"
                    << ", const labelList&, labelHashSet*): "
                    << "face " << face0 << " points the wrong way. " << endl
                    << "Pyramid volume: " << -pyrVolOwn
                    << " Face " << f[face0] << " area: " << f[face0].mag(p)
                    << " Owner cell: " << own[face0] << endl
                    << "Owner cell vertex labels: "
                    << mesh.cells()[own[face0]].labels(f)
                    << endl;
            }


            if (setPtr)
            {
                setPtr->insert(face0);
            }

            nErrorPyrs++;
        }

        // Create the neighbour pyramid - it will have positive volume
        scalar pyrVolNbr =
            pyramidPointFaceRef(f[face0], cellCentres[own[face1]]).mag(p);

        if (pyrVolNbr < minPyrVol)
        {
            if (report)
            {
                Pout<< "bool polyMeshGeometry::checkFacePyramids("
                    << "const bool, const scalar, const pointField&"
                    << ", const labelList&, labelHashSet*): "
                    << "face " << face0 << " points the wrong way. " << endl
                    << "Pyramid volume: " << -pyrVolNbr
                    << " Face " << f[face0] << " area: " << f[face0].mag(p)
                    << " Neighbour cell: " << own[face1] << endl
                    << "Neighbour cell vertex labels: "
                    << mesh.cells()[own[face1]].labels(f)
                    << endl;
            }

            if (setPtr)
            {
                setPtr->insert(face0);
            }

            nErrorPyrs++;
        }
    }

    reduce(nErrorPyrs, sumOp<label>());

    if (nErrorPyrs > 0)
    {
        if (report)
        {
            SeriousErrorInFunction
                << "Error in face pyramids: faces pointing the wrong way."
                << endl;
        }

        return true;
    }
    else
    {
        if (report)
        {
            Info<< "Face pyramids OK.\n" << endl;
        }

        return false;
    }
}


bool Foam::polyMeshGeometry::checkFaceTets
(
    const bool report,
    const scalar minTetQuality,
    const polyMesh& mesh,
    const vectorField& cellCentres,
    const vectorField& faceCentres,
    const pointField& p,
    const labelList& checkFaces,
    const List<labelPair>& baffles,
    labelHashSet* setPtr
)
{
    // check whether decomposing each cell into tets results in
    // positive volume, non-flat tets
    const labelList& own = mesh.faceOwner();
    const labelList& nei = mesh.faceNeighbour();
    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    // Calculate coupled cell centre
    pointField neiCc(mesh.nFaces() - mesh.nInternalFaces());

    for (label facei = mesh.nInternalFaces(); facei < mesh.nFaces(); facei++)
    {
        neiCc[facei - mesh.nInternalFaces()] = cellCentres[own[facei]];
    }

    syncTools::swapBoundaryFacePositions(mesh, neiCc);

    label nErrorTets = 0;

    forAll(checkFaces, i)
    {
        label facei = checkFaces[i];

        // Create the owner pyramid - note: exchange cell and face centre
        // to get positive volume.
        bool tetError = checkFaceTet
        (
            mesh,
            report,
            minTetQuality,
            p,
            facei,
            cellCentres[own[facei]],    // face centre
            faceCentres[facei],         // cell centre
            setPtr
        );

        if (tetError)
        {
            nErrorTets++;
        }

        if (mesh.isInternalFace(facei))
        {
            // Create the neighbour tets - they will have positive volume
            bool tetError = checkFaceTet
            (
                mesh,
                report,
                minTetQuality,
                p,
                facei,
                faceCentres[facei],         // face centre
                cellCentres[nei[facei]],    // cell centre
                setPtr
            );

            if (tetError)
            {
                nErrorTets++;
            }

            if
            (
                polyMeshTetDecomposition::findSharedBasePoint
                (
                    mesh,
                    facei,
                    minTetQuality,
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
            label patchi = patches.whichPatch(facei);

            if (patches[patchi].coupled())
            {
                if
                (
                    polyMeshTetDecomposition::findSharedBasePoint
                    (
                        mesh,
                        facei,
                        neiCc[facei - mesh.nInternalFaces()],
                        minTetQuality,
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
                if
                (
                    polyMeshTetDecomposition::findBasePoint
                    (
                        mesh,
                        facei,
                        minTetQuality,
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
        }
    }

    forAll(baffles, i)
    {
        label face0 = baffles[i].first();
        label face1 = baffles[i].second();

        bool tetError = checkFaceTet
        (
            mesh,
            report,
            minTetQuality,
            p,
            face0,
            cellCentres[own[face0]],    // face centre
            faceCentres[face0],         // cell centre
            setPtr
        );

        if (tetError)
        {
            nErrorTets++;
        }

        // Create the neighbour tets - they will have positive volume
        tetError = checkFaceTet
        (
            mesh,
            report,
            minTetQuality,
            p,
            face0,
            faceCentres[face0],         // face centre
            cellCentres[own[face1]],    // cell centre
            setPtr
        );

        if (tetError)
        {
            nErrorTets++;
        }

        if
        (
            polyMeshTetDecomposition::findSharedBasePoint
            (
                mesh,
                face0,
                cellCentres[own[face1]],
                minTetQuality,
                report
            ) == -1
        )
        {
            if (setPtr)
            {
                setPtr->insert(face0);
            }

            nErrorTets++;
        }
    }

    reduce(nErrorTets, sumOp<label>());

    if (nErrorTets > 0)
    {
        if (report)
        {
            SeriousErrorInFunction
                << "Error in face decomposition: negative tets."
                << endl;
        }

        return true;
    }
    else
    {
        if (report)
        {
            Info<< "Face tets OK.\n" << endl;
        }

        return false;
    }
}


bool Foam::polyMeshGeometry::checkFaceSkewness
(
    const bool report,
    const scalar internalSkew,
    const scalar boundarySkew,
    const polyMesh& mesh,
    const pointField& points,
    const vectorField& cellCentres,
    const vectorField& faceCentres,
    const vectorField& faceAreas,
    const labelList& checkFaces,
    const List<labelPair>& baffles,
    labelHashSet* setPtr
)
{
    // Warn if the skew correction vector is more than skew times
    // larger than the face area vector

    const labelList& own = mesh.faceOwner();
    const labelList& nei = mesh.faceNeighbour();
    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    // Calculate coupled cell centre
    pointField neiCc;
    syncTools::swapBoundaryCellPositions(mesh, cellCentres, neiCc);

    scalar maxSkew = 0;

    label nWarnSkew = 0;

    forAll(checkFaces, i)
    {
        label facei = checkFaces[i];

        if (mesh.isInternalFace(facei))
        {
            scalar skewness = primitiveMeshTools::faceSkewness
            (
                mesh,
                points,
                faceCentres,
                faceAreas,

                facei,
                cellCentres[own[facei]],
                cellCentres[nei[facei]]
            );

            // Check if the skewness vector is greater than the PN vector.
            // This does not cause trouble but is a good indication of a poor
            // mesh.
            if (skewness > internalSkew)
            {
                if (report)
                {
                    Pout<< "Severe skewness for face " << facei
                        << " skewness = " << skewness << endl;
                }

                if (setPtr)
                {
                    setPtr->insert(facei);
                }

                nWarnSkew++;
            }

            maxSkew = max(maxSkew, skewness);
        }
        else if (patches[patches.whichPatch(facei)].coupled())
        {
            scalar skewness = primitiveMeshTools::faceSkewness
            (
                mesh,
                points,
                faceCentres,
                faceAreas,

                facei,
                cellCentres[own[facei]],
                neiCc[facei-mesh.nInternalFaces()]
            );

            // Check if the skewness vector is greater than the PN vector.
            // This does not cause trouble but is a good indication of a poor
            // mesh.
            if (skewness > internalSkew)
            {
                if (report)
                {
                    Pout<< "Severe skewness for coupled face " << facei
                        << " skewness = " << skewness << endl;
                }

                if (setPtr)
                {
                    setPtr->insert(facei);
                }

                nWarnSkew++;
            }

            maxSkew = max(maxSkew, skewness);
        }
        else
        {
            scalar skewness = primitiveMeshTools::boundaryFaceSkewness
            (
                mesh,
                points,
                faceCentres,
                faceAreas,

                facei,
                cellCentres[own[facei]]
            );


            // Check if the skewness vector is greater than the PN vector.
            // This does not cause trouble but is a good indication of a poor
            // mesh.
            if (skewness > boundarySkew)
            {
                if (report)
                {
                    Pout<< "Severe skewness for boundary face " << facei
                        << " skewness = " << skewness << endl;
                }

                if (setPtr)
                {
                    setPtr->insert(facei);
                }

                nWarnSkew++;
            }

            maxSkew = max(maxSkew, skewness);
        }
    }

    forAll(baffles, i)
    {
        label face0 = baffles[i].first();
        label face1 = baffles[i].second();

        const point& ownCc = cellCentres[own[face0]];
        const point& neiCc = cellCentres[own[face1]];

        scalar skewness = primitiveMeshTools::faceSkewness
        (
            mesh,
            points,
            faceCentres,
            faceAreas,

            face0,
            ownCc,
            neiCc
        );

        // Check if the skewness vector is greater than the PN vector.
        // This does not cause trouble but is a good indication of a poor
        // mesh.
        if (skewness > internalSkew)
        {
            if (report)
            {
                Pout<< "Severe skewness for face " << face0
                    << " skewness = " << skewness << endl;
            }

            if (setPtr)
            {
                setPtr->insert(face0);
            }

            nWarnSkew++;
        }

        maxSkew = max(maxSkew, skewness);
    }


    reduce(maxSkew, maxOp<scalar>());
    reduce(nWarnSkew, sumOp<label>());

    if (nWarnSkew > 0)
    {
        if (report)
        {
            WarningInFunction
                << 100*maxSkew
                << " percent.\nThis may impair the quality of the result." << nl
                << nWarnSkew << " highly skew faces detected."
                << endl;
        }

        return true;
    }
    else
    {
        if (report)
        {
            Info<< "Max skewness = " << 100*maxSkew
                << " percent.  Face skewness OK.\n" << endl;
        }

        return false;
    }
}


bool Foam::polyMeshGeometry::checkFaceWeights
(
    const bool report,
    const scalar warnWeight,
    const polyMesh& mesh,
    const vectorField& cellCentres,
    const vectorField& faceCentres,
    const vectorField& faceAreas,
    const labelList& checkFaces,
    const List<labelPair>& baffles,
    labelHashSet* setPtr
)
{
    // Warn if the delta factor (0..1) is too large.

    const labelList& own = mesh.faceOwner();
    const labelList& nei = mesh.faceNeighbour();
    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    // Calculate coupled cell centre
    pointField neiCc(mesh.nFaces()-mesh.nInternalFaces());

    for (label facei = mesh.nInternalFaces(); facei < mesh.nFaces(); facei++)
    {
        neiCc[facei-mesh.nInternalFaces()] = cellCentres[own[facei]];
    }
    syncTools::swapBoundaryFacePositions(mesh, neiCc);


    scalar minWeight = great;

    label nWarnWeight = 0;

    forAll(checkFaces, i)
    {
        label facei = checkFaces[i];

        const point& fc = faceCentres[facei];
        const vector& fa = faceAreas[facei];

        scalar dOwn = mag(fa & (fc-cellCentres[own[facei]]));

        if (mesh.isInternalFace(facei))
        {
            scalar dNei = mag(fa & (cellCentres[nei[facei]]-fc));
            scalar weight = min(dNei,dOwn)/(dNei+dOwn+vSmall);

            if (weight < warnWeight)
            {
                if (report)
                {
                    Pout<< "Small weighting factor for face " << facei
                        << " weight = " << weight << endl;
                }

                if (setPtr)
                {
                    setPtr->insert(facei);
                }

                nWarnWeight++;
            }

            minWeight = min(minWeight, weight);
        }
        else
        {
            label patchi = patches.whichPatch(facei);

            if (patches[patchi].coupled())
            {
                scalar dNei = mag(fa & (neiCc[facei-mesh.nInternalFaces()]-fc));
                scalar weight = min(dNei,dOwn)/(dNei+dOwn+vSmall);

                if (weight < warnWeight)
                {
                    if (report)
                    {
                        Pout<< "Small weighting factor for face " << facei
                            << " weight = " << weight << endl;
                    }

                    if (setPtr)
                    {
                        setPtr->insert(facei);
                    }

                    nWarnWeight++;
                }

                minWeight = min(minWeight, weight);
            }
        }
    }

    forAll(baffles, i)
    {
        label face0 = baffles[i].first();
        label face1 = baffles[i].second();

        const point& ownCc = cellCentres[own[face0]];
        const point& fc = faceCentres[face0];
        const vector& fa = faceAreas[face0];

        scalar dOwn = mag(fa & (fc-ownCc));
        scalar dNei = mag(fa & (cellCentres[own[face1]]-fc));
        scalar weight = min(dNei,dOwn)/(dNei+dOwn+vSmall);

        if (weight < warnWeight)
        {
            if (report)
            {
                Pout<< "Small weighting factor for face " << face0
                    << " weight = " << weight << endl;
            }

            if (setPtr)
            {
                setPtr->insert(face0);
            }

            nWarnWeight++;
        }

        minWeight = min(minWeight, weight);
    }

    reduce(minWeight, minOp<scalar>());
    reduce(nWarnWeight, sumOp<label>());

    if (minWeight < warnWeight)
    {
        if (report)
        {
            WarningInFunction
                << minWeight << '.' << nl
                << nWarnWeight << " faces with small weights detected."
                << endl;
        }

        return true;
    }
    else
    {
        if (report)
        {
            Info<< "Min weight = " << minWeight
                << ".  Weights OK.\n" << endl;
        }

        return false;
    }
}


bool Foam::polyMeshGeometry::checkVolRatio
(
    const bool report,
    const scalar warnRatio,
    const polyMesh& mesh,
    const scalarField& cellVolumes,
    const labelList& checkFaces,
    const List<labelPair>& baffles,
    labelHashSet* setPtr
)
{
    // Warn if the volume ratio between neighbouring cells is too large

    const labelList& own = mesh.faceOwner();
    const labelList& nei = mesh.faceNeighbour();
    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    // Calculate coupled cell vol
    scalarField neiVols(mesh.nFaces()-mesh.nInternalFaces());

    for (label facei = mesh.nInternalFaces(); facei < mesh.nFaces(); facei++)
    {
        neiVols[facei-mesh.nInternalFaces()] = cellVolumes[own[facei]];
    }
    syncTools::swapBoundaryFaceList(mesh, neiVols);


    scalar minRatio = great;

    label nWarnRatio = 0;

    forAll(checkFaces, i)
    {
        label facei = checkFaces[i];

        scalar ownVol = mag(cellVolumes[own[facei]]);

        scalar neiVol = -great;

        if (mesh.isInternalFace(facei))
        {
            neiVol = mag(cellVolumes[nei[facei]]);
        }
        else
        {
            label patchi = patches.whichPatch(facei);

            if (patches[patchi].coupled())
            {
                neiVol = mag(neiVols[facei-mesh.nInternalFaces()]);
            }
        }

        if (neiVol >= 0)
        {
            scalar ratio = min(ownVol, neiVol) / (max(ownVol, neiVol) + vSmall);

            if (ratio < warnRatio)
            {
                if (report)
                {
                    Pout<< "Small ratio for face " << facei
                        << " ratio = " << ratio << endl;
                }

                if (setPtr)
                {
                    setPtr->insert(facei);
                }

                nWarnRatio++;
            }

            minRatio = min(minRatio, ratio);
        }
    }

    forAll(baffles, i)
    {
        label face0 = baffles[i].first();
        label face1 = baffles[i].second();

        scalar ownVol = mag(cellVolumes[own[face0]]);

        scalar neiVol = mag(cellVolumes[own[face1]]);

        if (neiVol >= 0)
        {
            scalar ratio = min(ownVol, neiVol) / (max(ownVol, neiVol) + vSmall);

            if (ratio < warnRatio)
            {
                if (report)
                {
                    Pout<< "Small ratio for face " << face0
                        << " ratio = " << ratio << endl;
                }

                if (setPtr)
                {
                    setPtr->insert(face0);
                }

                nWarnRatio++;
            }

            minRatio = min(minRatio, ratio);
        }
    }

    reduce(minRatio, minOp<scalar>());
    reduce(nWarnRatio, sumOp<label>());

    if (minRatio < warnRatio)
    {
        if (report)
        {
            WarningInFunction
                << minRatio << '.' << nl
                << nWarnRatio << " faces with small ratios detected."
                << endl;
        }

        return true;
    }
    else
    {
        if (report)
        {
            Info<< "Min ratio = " << minRatio
                << ".  Ratios OK.\n" << endl;
        }

        return false;
    }
}


// Check convexity of angles in a face. Allow a slight non-convexity.
// E.g. maxDeg = 10 allows for angles < 190 (or 10 degrees concavity)
// (if truly concave and points not visible from face centre the face-pyramid
//  check in checkMesh will fail)
bool Foam::polyMeshGeometry::checkFaceAngles
(
    const bool report,
    const scalar maxDeg,
    const polyMesh& mesh,
    const vectorField& faceAreas,
    const pointField& p,
    const labelList& checkFaces,
    labelHashSet* setPtr
)
{
    if (maxDeg < -small || maxDeg > 180+small)
    {
        FatalErrorInFunction
            << "maxDeg should be [0..180] but is now " << maxDeg
            << abort(FatalError);
    }

    const scalar maxSin = Foam::sin(degToRad(maxDeg));

    const faceList& fcs = mesh.faces();

    scalar maxEdgeSin = 0.0;

    label nConcave = 0;

    label errorFacei = -1;

    forAll(checkFaces, i)
    {
        label facei = checkFaces[i];

        const face& f = fcs[facei];

        vector faceNormal = faceAreas[facei];
        faceNormal /= mag(faceNormal) + vSmall;

        // Get edge from f[0] to f[size-1];
        vector ePrev(p[f.first()] - p[f.last()]);
        scalar magEPrev = mag(ePrev);
        ePrev /= magEPrev + vSmall;

        forAll(f, fp0)
        {
            // Get vertex after fp
            label fp1 = f.fcIndex(fp0);

            // Normalized vector between two consecutive points
            vector e10(p[f[fp1]] - p[f[fp0]]);
            scalar magE10 = mag(e10);
            e10 /= magE10 + vSmall;

            if (magEPrev > small && magE10 > small)
            {
                vector edgeNormal = ePrev ^ e10;
                scalar magEdgeNormal = mag(edgeNormal);

                if (magEdgeNormal < maxSin)
                {
                    // Edges (almost) aligned -> face is ok.
                }
                else
                {
                    // Check normal
                    edgeNormal /= magEdgeNormal;

                    if ((edgeNormal & faceNormal) < small)
                    {
                        if (facei != errorFacei)
                        {
                            // Count only one error per face.
                            errorFacei = facei;
                            nConcave++;
                        }

                        if (setPtr)
                        {
                            setPtr->insert(facei);
                        }

                        maxEdgeSin = max(maxEdgeSin, magEdgeNormal);
                    }
                }
            }

            ePrev = e10;
            magEPrev = magE10;
        }
    }

    reduce(nConcave, sumOp<label>());
    reduce(maxEdgeSin, maxOp<scalar>());

    if (report)
    {
        if (maxEdgeSin > small)
        {
            scalar maxConcaveDegr =
                radToDeg(Foam::asin(Foam::min(1.0, maxEdgeSin)));

            Info<< "There are " << nConcave
                << " faces with concave angles between consecutive"
                << " edges. Max concave angle = " << maxConcaveDegr
                << " degrees.\n" << endl;
        }
        else
        {
            Info<< "All angles in faces are convex or less than "  << maxDeg
                << " degrees concave.\n" << endl;
        }
    }

    if (nConcave > 0)
    {
        if (report)
        {
            WarningInFunction
                << nConcave  << " face points with severe concave angle (> "
                << maxDeg << " deg) found.\n"
                << endl;
        }

        return true;
    }
    else
    {
        return false;
    }
}


// Check twist of faces. Is calculated as the difference between normals of
// individual triangles and the cell-cell centre edge.
bool Foam::polyMeshGeometry::checkFaceTwist
(
    const bool report,
    const scalar minTwist,
    const polyMesh& mesh,
    const vectorField& cellCentres,
    const vectorField& faceAreas,
    const vectorField& faceCentres,
    const pointField& p,
    const labelList& checkFaces,
    labelHashSet* setPtr
)
{
    if (minTwist < -1-small || minTwist > 1+small)
    {
        FatalErrorInFunction
            << "minTwist should be [-1..1] but is now " << minTwist
            << abort(FatalError);
    }


    const faceList& fcs = mesh.faces();


    label nWarped = 0;

//    forAll(checkFaces, i)
//    {
//        label facei = checkFaces[i];
//
//        const face& f = fcs[facei];
//
//        scalar magArea = mag(faceAreas[facei]);
//
//        if (f.size() > 3 && magArea > vSmall)
//        {
//            const vector nf = faceAreas[facei] / magArea;
//
//            const point& fc = faceCentres[facei];
//
//            forAll(f, fpI)
//            {
//                vector triArea
//                (
//                    triPointRef
//                    (
//                        p[f[fpI]],
//                        p[f.nextLabel(fpI)],
//                        fc
//                    ).area()
//                );
//
//                scalar magTri = mag(triArea);
//
//                if (magTri > vSmall && ((nf & triArea/magTri) < minTwist))
//                {
//                    nWarped++;
//
//                    if (setPtr)
//                    {
//                        setPtr->insert(facei);
//                    }
//
//                    break;
//                }
//            }
//        }
//    }


    const labelList& own = mesh.faceOwner();
    const labelList& nei = mesh.faceNeighbour();
    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    // Calculate coupled cell centre
    pointField neiCc(mesh.nFaces()-mesh.nInternalFaces());

    for (label facei = mesh.nInternalFaces(); facei < mesh.nFaces(); facei++)
    {
        neiCc[facei-mesh.nInternalFaces()] = cellCentres[own[facei]];
    }
    syncTools::swapBoundaryFacePositions(mesh, neiCc);

    forAll(checkFaces, i)
    {
        label facei = checkFaces[i];

        const face& f = fcs[facei];

        if (f.size() > 3)
        {
            vector nf(Zero);

            if (mesh.isInternalFace(facei))
            {
                nf = cellCentres[nei[facei]] - cellCentres[own[facei]];
                nf /= mag(nf) + vSmall;
            }
            else if (patches[patches.whichPatch(facei)].coupled())
            {
                nf =
                    neiCc[facei-mesh.nInternalFaces()]
                  - cellCentres[own[facei]];
                nf /= mag(nf) + vSmall;
            }
            else
            {
                nf = faceCentres[facei] - cellCentres[own[facei]];
                nf /= mag(nf) + vSmall;
            }

            if (nf != vector::zero)
            {
                const point& fc = faceCentres[facei];

                forAll(f, fpI)
                {
                    vector triArea
                    (
                        triPointRef
                        (
                            p[f[fpI]],
                            p[f.nextLabel(fpI)],
                            fc
                        ).area()
                    );

                    scalar magTri = mag(triArea);

                    if (magTri > vSmall && ((nf & triArea/magTri) < minTwist))
                    {
                        nWarped++;

                        if (setPtr)
                        {
                            setPtr->insert(facei);
                        }

                        break;
                    }
                }
            }
        }
    }

    reduce(nWarped, sumOp<label>());

    if (report)
    {
        if (nWarped> 0)
        {
            Info<< "There are " << nWarped
                << " faces with cosine of the angle"
                << " between triangle normal and face normal less than "
                << minTwist << nl << endl;
        }
        else
        {
            Info<< "All faces are flat in that the cosine of the angle"
                << " between triangle normal and face normal less than "
                << minTwist << nl << endl;
        }
    }

    if (nWarped > 0)
    {
        if (report)
        {
            WarningInFunction
                << nWarped  << " faces with severe warpage "
                << "(cosine of the angle between triangle normal and "
                << "face normal < " << minTwist << ") found.\n"
                << endl;
        }

        return true;
    }
    else
    {
        return false;
    }
}


// Like checkFaceTwist but compares normals of consecutive triangles.
bool Foam::polyMeshGeometry::checkTriangleTwist
(
    const bool report,
    const scalar minTwist,
    const polyMesh& mesh,
    const vectorField& faceAreas,
    const vectorField& faceCentres,
    const pointField& p,
    const labelList& checkFaces,
    labelHashSet* setPtr
)
{
    if (minTwist < -1-small || minTwist > 1+small)
    {
        FatalErrorInFunction
            << "minTwist should be [-1..1] but is now " << minTwist
            << abort(FatalError);
    }

    const faceList& fcs = mesh.faces();

    label nWarped = 0;

    forAll(checkFaces, i)
    {
        label facei = checkFaces[i];

        const face& f = fcs[facei];

        if (f.size() > 3)
        {
            const point& fc = faceCentres[facei];

            // Find starting triangle (at startFp) with non-zero area
            label startFp = -1;
            vector prevN;

            forAll(f, fp)
            {
                prevN = triPointRef
                (
                    p[f[fp]],
                    p[f.nextLabel(fp)],
                    fc
                ).area();

                scalar magTri = mag(prevN);

                if (magTri > vSmall)
                {
                    startFp = fp;
                    prevN /= magTri;
                    break;
                }
            }

            if (startFp != -1)
            {
                label fp = startFp;

                do
                {
                    fp = f.fcIndex(fp);

                    vector triN
                    (
                        triPointRef
                        (
                            p[f[fp]],
                            p[f.nextLabel(fp)],
                            fc
                        ).area()
                    );
                    scalar magTri = mag(triN);

                    if (magTri > vSmall)
                    {
                        triN /= magTri;

                        if ((prevN & triN) < minTwist)
                        {
                            nWarped++;

                            if (setPtr)
                            {
                                setPtr->insert(facei);
                            }

                            break;
                        }

                        prevN = triN;
                    }
                    else if (minTwist > 0)
                    {
                        nWarped++;

                        if (setPtr)
                        {
                            setPtr->insert(facei);
                        }

                        break;
                    }
                }
                while (fp != startFp);
            }
        }
    }


    reduce(nWarped, sumOp<label>());

    if (report)
    {
        if (nWarped> 0)
        {
            Info<< "There are " << nWarped
                << " faces with cosine of the angle"
                << " between consecutive triangle normals less than "
                << minTwist << nl << endl;
        }
        else
        {
            Info<< "All faces are flat in that the cosine of the angle"
                << " between consecutive triangle normals is less than "
                << minTwist << nl << endl;
        }
    }

    if (nWarped > 0)
    {
        if (report)
        {
            WarningInFunction
                << nWarped  << " faces with severe warpage "
                << "(cosine of the angle between consecutive triangle normals"
                << " < " << minTwist << ") found.\n"
                << endl;
        }

        return true;
    }
    else
    {
        return false;
    }
}


bool Foam::polyMeshGeometry::checkFaceFlatness
(
    const bool report,
    const scalar minFlatness,
    const polyMesh& mesh,
    const vectorField& faceAreas,
    const vectorField& faceCentres,
    const pointField& p,
    const labelList& checkFaces,
    labelHashSet* setPtr
)
{
    if (minFlatness < -small || minFlatness > 1+small)
    {
        FatalErrorInFunction
            << "minFlatness should be [0..1] but is now " << minFlatness
            << abort(FatalError);
    }

    const faceList& fcs = mesh.faces();

    label nWarped = 0;

    forAll(checkFaces, i)
    {
        label facei = checkFaces[i];

        const face& f = fcs[facei];

        if (f.size() > 3)
        {
            const point& fc = faceCentres[facei];

            // Sum triangle areas
            scalar sumArea = 0.0;

            forAll(f, fp)
            {
                sumArea += triPointRef
                (
                    p[f[fp]],
                    p[f.nextLabel(fp)],
                    fc
                ).mag();
            }

            if (sumArea/mag(faceAreas[facei]) < minFlatness)
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

    if (report)
    {
        if (nWarped> 0)
        {
            Info<< "There are " << nWarped
                << " faces with area of individual triangles"
                << " compared to overall area less than "
                << minFlatness << nl << endl;
        }
        else
        {
            Info<< "All faces are flat in that the area of individual triangles"
                << " compared to overall area is less than "
                << minFlatness << nl << endl;
        }
    }

    if (nWarped > 0)
    {
        if (report)
        {
            WarningInFunction
                << nWarped  << " non-flat faces "
                << "(area of individual triangles"
                << " compared to overall area"
                << " < " << minFlatness << ") found.\n"
                << endl;
        }

        return true;
    }
    else
    {
        return false;
    }
}


bool Foam::polyMeshGeometry::checkFaceArea
(
    const bool report,
    const scalar minArea,
    const polyMesh& mesh,
    const vectorField& faceAreas,
    const labelList& checkFaces,
    labelHashSet* setPtr
)
{
    label nZeroArea = 0;

    forAll(checkFaces, i)
    {
        label facei = checkFaces[i];

        if (mag(faceAreas[facei]) < minArea)
        {
            if (setPtr)
            {
                setPtr->insert(facei);
            }
            nZeroArea++;
        }
    }


    reduce(nZeroArea, sumOp<label>());

    if (report)
    {
        if (nZeroArea > 0)
        {
            Info<< "There are " << nZeroArea
                << " faces with area < " << minArea << '.' << nl << endl;
        }
        else
        {
            Info<< "All faces have area > " << minArea << '.' << nl << endl;
        }
    }

    if (nZeroArea > 0)
    {
        if (report)
        {
            WarningInFunction
                << nZeroArea  << " faces with area < " << minArea
                << " found.\n"
                << endl;
        }

        return true;
    }
    else
    {
        return false;
    }
}


bool Foam::polyMeshGeometry::checkCellDeterminant
(
    const bool report,
    const scalar warnDet,
    const polyMesh& mesh,
    const vectorField& faceAreas,
    const labelList& checkFaces,
    const labelList& affectedCells,
    labelHashSet* setPtr
)
{
    const cellList& cells = mesh.cells();

    scalar minDet = great;
    scalar sumDet = 0.0;
    label nSumDet = 0;
    label nWarnDet = 0;

    forAll(affectedCells, i)
    {
        const cell& cFaces = cells[affectedCells[i]];

        tensor areaSum(Zero);
        scalar magAreaSum = 0;

        forAll(cFaces, cFacei)
        {
            label facei = cFaces[cFacei];

            scalar magArea = mag(faceAreas[facei]);

            magAreaSum += magArea;
            areaSum += faceAreas[facei]*(faceAreas[facei]/(magArea+vSmall));
        }

        scalar scaledDet = det(areaSum/(magAreaSum+vSmall))/0.037037037037037;

        minDet = min(minDet, scaledDet);
        sumDet += scaledDet;
        nSumDet++;

        if (scaledDet < warnDet)
        {
            if (setPtr)
            {
                // Insert all faces of the cell.
                forAll(cFaces, cFacei)
                {
                    label facei = cFaces[cFacei];
                    setPtr->insert(facei);
                }
            }
            nWarnDet++;
        }
    }

    reduce(minDet, minOp<scalar>());
    reduce(sumDet, sumOp<scalar>());
    reduce(nSumDet, sumOp<label>());
    reduce(nWarnDet, sumOp<label>());

    if (report)
    {
        if (nSumDet > 0)
        {
            Info<< "Cell determinant (1 = uniform cube) : average = "
                << sumDet / nSumDet << "  min = " << minDet << endl;
        }

        if (nWarnDet > 0)
        {
            Info<< "There are " << nWarnDet
                << " cells with determinant < " << warnDet << '.' << nl
                << endl;
        }
        else
        {
            Info<< "All faces have determinant > " << warnDet << '.' << nl
                << endl;
        }
    }

    if (nWarnDet > 0)
    {
        if (report)
        {
            WarningInFunction
                << nWarnDet << " cells with determinant < " << warnDet
                << " found.\n"
                << endl;
        }

        return true;
    }
    else
    {
        return false;
    }
}


bool Foam::polyMeshGeometry::checkFaceDotProduct
(
    const bool report,
    const scalar orthWarn,
    const labelList& checkFaces,
    const List<labelPair>& baffles,
    labelHashSet* setPtr
) const
{
    return checkFaceDotProduct
    (
        report,
        orthWarn,
        mesh_,
        cellCentres_,
        faceAreas_,
        checkFaces,
        baffles,
        setPtr
    );
}


bool Foam::polyMeshGeometry::checkFacePyramids
(
    const bool report,
    const scalar minPyrVol,
    const pointField& p,
    const labelList& checkFaces,
    const List<labelPair>& baffles,
    labelHashSet* setPtr
) const
{
    return checkFacePyramids
    (
        report,
        minPyrVol,
        mesh_,
        cellCentres_,
        p,
        checkFaces,
        baffles,
        setPtr
    );
}


bool Foam::polyMeshGeometry::checkFaceTets
(
    const bool report,
    const scalar minTetQuality,
    const pointField& p,
    const labelList& checkFaces,
    const List<labelPair>& baffles,
    labelHashSet* setPtr
) const
{
    return checkFaceTets
    (
        report,
        minTetQuality,
        mesh_,
        cellCentres_,
        faceCentres_,
        p,
        checkFaces,
        baffles,
        setPtr
    );
}


//bool Foam::polyMeshGeometry::checkFaceSkewness
//(
//    const bool report,
//    const scalar internalSkew,
//    const scalar boundarySkew,
//    const labelList& checkFaces,
//    const List<labelPair>& baffles,
//    labelHashSet* setPtr
//) const
//{
//    return checkFaceSkewness
//    (
//        report,
//        internalSkew,
//        boundarySkew,
//        mesh_,
//        cellCentres_,
//        faceCentres_,
//        faceAreas_,
//        checkFaces,
//        baffles,
//        setPtr
//    );
//}


bool Foam::polyMeshGeometry::checkFaceWeights
(
    const bool report,
    const scalar warnWeight,
    const labelList& checkFaces,
    const List<labelPair>& baffles,
    labelHashSet* setPtr
) const
{
    return checkFaceWeights
    (
        report,
        warnWeight,
        mesh_,
        cellCentres_,
        faceCentres_,
        faceAreas_,
        checkFaces,
        baffles,
        setPtr
    );
}


bool Foam::polyMeshGeometry::checkVolRatio
(
    const bool report,
    const scalar warnRatio,
    const labelList& checkFaces,
    const List<labelPair>& baffles,
    labelHashSet* setPtr
) const
{
    return checkVolRatio
    (
        report,
        warnRatio,
        mesh_,
        cellVolumes_,
        checkFaces,
        baffles,
        setPtr
    );
}


bool Foam::polyMeshGeometry::checkFaceAngles
(
    const bool report,
    const scalar maxDeg,
    const pointField& p,
    const labelList& checkFaces,
    labelHashSet* setPtr
) const
{
    return checkFaceAngles
    (
        report,
        maxDeg,
        mesh_,
        faceAreas_,
        p,
        checkFaces,
        setPtr
    );
}


bool Foam::polyMeshGeometry::checkFaceTwist
(
    const bool report,
    const scalar minTwist,
    const pointField& p,
    const labelList& checkFaces,
    labelHashSet* setPtr
) const
{
    return checkFaceTwist
    (
        report,
        minTwist,
        mesh_,
        cellCentres_,
        faceAreas_,
        faceCentres_,
        p,
        checkFaces,
        setPtr
    );
}


bool Foam::polyMeshGeometry::checkTriangleTwist
(
    const bool report,
    const scalar minTwist,
    const pointField& p,
    const labelList& checkFaces,
    labelHashSet* setPtr
) const
{
    return checkTriangleTwist
    (
        report,
        minTwist,
        mesh_,
        faceAreas_,
        faceCentres_,
        p,
        checkFaces,
        setPtr
    );
}


bool Foam::polyMeshGeometry::checkFaceFlatness
(
    const bool report,
    const scalar minFlatness,
    const pointField& p,
    const labelList& checkFaces,
    labelHashSet* setPtr
) const
{
    return checkFaceFlatness
    (
        report,
        minFlatness,
        mesh_,
        faceAreas_,
        faceCentres_,
        p,
        checkFaces,
        setPtr
    );
}


bool Foam::polyMeshGeometry::checkFaceArea
(
    const bool report,
    const scalar minArea,
    const labelList& checkFaces,
    labelHashSet* setPtr
) const
{
    return checkFaceArea
    (
        report,
        minArea,
        mesh_,
        faceAreas_,
        checkFaces,
        setPtr
    );
}


bool Foam::polyMeshGeometry::checkCellDeterminant
(
    const bool report,
    const scalar warnDet,
    const labelList& checkFaces,
    const labelList& affectedCells,
    labelHashSet* setPtr
) const
{
    return checkCellDeterminant
    (
        report,
        warnDet,
        mesh_,
        faceAreas_,
        checkFaces,
        affectedCells,
        setPtr
    );
}


// ************************************************************************* //
