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

#include "polyMeshGeometry.H"
#include "polyMeshTetDecomposition.H"
#include "pyramidPointFaceRef.H"
#include "tetrahedron.H"
#include "syncTools.H"
#include "unitConversion.H"
#include "primitiveMeshTools.H"

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

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
            vector sumN = vector::zero;
            scalar sumA = 0.0;
            vector sumAc = vector::zero;

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

            faceCentres_[facei] = (1.0/3.0)*sumAc/(sumA + VSMALL);
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
    UIndirectList<vector>(cellCentres_, changedCells) = vector::zero;
    UIndirectList<scalar>(cellVolumes_, changedCells) = 0.0;

    const labelList& own = mesh_.faceOwner();
    const labelList& nei = mesh_.faceNeighbour();

    // first estimate the approximate cell centre as the average of face centres

    vectorField cEst(mesh_.nCells());
    UIndirectList<vector>(cEst, changedCells) = vector::zero;
    scalarField nCellFaces(mesh_.nCells());
    UIndirectList<scalar>(nCellFaces, changedCells) = 0.0;

    forAll(changedFaces, i)
    {
        label faceI = changedFaces[i];
        cEst[own[faceI]] += faceCentres_[faceI];
        nCellFaces[own[faceI]] += 1;

        if (mesh_.isInternalFace(faceI))
        {
            cEst[nei[faceI]] += faceCentres_[faceI];
            nCellFaces[nei[faceI]] += 1;
        }
    }

    forAll(changedCells, i)
    {
        label cellI = changedCells[i];
        cEst[cellI] /= nCellFaces[cellI];
    }

    forAll(changedFaces, i)
    {
        label faceI = changedFaces[i];

        // Calculate 3*face-pyramid volume
        scalar pyr3Vol = max
        (
            faceAreas_[faceI] & (faceCentres_[faceI] - cEst[own[faceI]]),
            VSMALL
        );

        // Calculate face-pyramid centre
        vector pc = (3.0/4.0)*faceCentres_[faceI] + (1.0/4.0)*cEst[own[faceI]];

        // Accumulate volume-weighted face-pyramid centre
        cellCentres_[own[faceI]] += pyr3Vol*pc;

        // Accumulate face-pyramid volume
        cellVolumes_[own[faceI]] += pyr3Vol;

        if (mesh_.isInternalFace(faceI))
        {
            // Calculate 3*face-pyramid volume
            scalar pyr3Vol = max
            (
                faceAreas_[faceI] & (cEst[nei[faceI]] - faceCentres_[faceI]),
                VSMALL
            );

            // Calculate face-pyramid centre
            vector pc =
                (3.0/4.0)*faceCentres_[faceI]
              + (1.0/4.0)*cEst[nei[faceI]];

            // Accumulate volume-weighted face-pyramid centre
            cellCentres_[nei[faceI]] += pyr3Vol*pc;

            // Accumulate face-pyramid volume
            cellVolumes_[nei[faceI]] += pyr3Vol;
        }
    }

    forAll(changedCells, i)
    {
        label cellI = changedCells[i];

        cellCentres_[cellI] /= cellVolumes_[cellI] + VSMALL;
        cellVolumes_[cellI] *= (1.0/3.0);
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
        label faceI = changedFaces[i];

        affectedCells.insert(own[faceI]);

        if (mesh.isInternalFace(faceI))
        {
            affectedCells.insert(nei[faceI]);
        }
    }
    return affectedCells.toc();
}


Foam::scalar Foam::polyMeshGeometry::checkNonOrtho
(
    const polyMesh& mesh,
    const bool report,
    const scalar severeNonorthogonalityThreshold,
    const label faceI,
    const vector& s,    // face area vector
    const vector& d,    // cc-cc vector

    label& severeNonOrth,
    label& errorNonOrth,
    labelHashSet* setPtr
)
{
    scalar dDotS = (d & s)/(mag(d)*mag(s) + VSMALL);

    if (dDotS < severeNonorthogonalityThreshold)
    {
        label nei = -1;

        if (mesh.isInternalFace(faceI))
        {
            nei = mesh.faceNeighbour()[faceI];
        }

        if (dDotS > SMALL)
        {
            if (report)
            {
                // Severe non-orthogonality but mesh still OK
                Pout<< "Severe non-orthogonality for face " << faceI
                    << " between cells " << mesh.faceOwner()[faceI]
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
                WarningIn
                (
                    "polyMeshGeometry::checkFaceDotProduct"
                    "(const bool, const scalar, const labelList&"
                    ", labelHashSet*)"
                )   << "Severe non-orthogonality detected for face "
                    << faceI
                    << " between cells " << mesh.faceOwner()[faceI]
                    << " and " << nei
                    << ": Angle = "
                    << radToDeg(::acos(dDotS))
                    << " deg." << endl;
            }

            errorNonOrth++;
        }

        if (setPtr)
        {
            setPtr->insert(faceI);
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
    const label faceI,
    const point& fc,    // face centre
    const point& cc,    // cell centre

    labelHashSet* setPtr
)
{
    const face& f = mesh.faces()[faceI];

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
                    << "face " << faceI
                    << " has a triangle that points the wrong way."
                     << endl
                    << "Tet quality: " << tetQual
                    << " Face " << faceI
                    << endl;
            }

            if (setPtr)
            {
                setPtr->insert(faceI);
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

//- Take over properties from mesh
void Foam::polyMeshGeometry::correct()
{
    faceAreas_ = mesh_.faceAreas();
    faceCentres_ = mesh_.faceCentres();
    cellCentres_ = mesh_.cellCentres();
    cellVolumes_ = mesh_.cellVolumes();
}


//- Recalculate on selected faces
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

    for (label faceI = mesh.nInternalFaces(); faceI < mesh.nFaces(); faceI++)
    {
        neiCc[faceI-mesh.nInternalFaces()] = cellCentres[own[faceI]];
    }

    syncTools::swapBoundaryFacePositions(mesh, neiCc);

    scalar minDDotS = GREAT;

    scalar sumDDotS = 0;
    label nDDotS = 0;

    label severeNonOrth = 0;

    label errorNonOrth = 0;

    forAll(checkFaces, i)
    {
        label faceI = checkFaces[i];

        const point& ownCc = cellCentres[own[faceI]];

        if (mesh.isInternalFace(faceI))
        {
            scalar dDotS = checkNonOrtho
            (
                mesh,
                report,
                severeNonorthogonalityThreshold,
                faceI,
                faceAreas[faceI],
                cellCentres[nei[faceI]] - ownCc,

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
            label patchI = patches.whichPatch(faceI);

            if (patches[patchI].coupled())
            {
                scalar dDotS = checkNonOrtho
                (
                    mesh,
                    report,
                    severeNonorthogonalityThreshold,
                    faceI,
                    faceAreas[faceI],
                    neiCc[faceI-mesh.nInternalFaces()] - ownCc,

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
            SeriousErrorIn
            (
                "polyMeshGeometry::checkFaceDotProduct"
                "(const bool, const scalar, const labelList&, labelHashSet*)"
            )   << "Error in non-orthogonality detected" << endl;
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
        label faceI = checkFaces[i];

        // Create the owner pyramid - it will have negative volume
        scalar pyrVol = pyramidPointFaceRef
        (
            f[faceI],
            cellCentres[own[faceI]]
        ).mag(p);

        if (pyrVol > -minPyrVol)
        {
            if (report)
            {
                Pout<< "bool polyMeshGeometry::checkFacePyramids("
                    << "const bool, const scalar, const pointField&"
                    << ", const labelList&, labelHashSet*): "
                    << "face " << faceI << " points the wrong way. " << endl
                    << "Pyramid volume: " << -pyrVol
                    << " Face " << f[faceI] << " area: " << f[faceI].mag(p)
                    << " Owner cell: " << own[faceI] << endl
                    << "Owner cell vertex labels: "
                    << mesh.cells()[own[faceI]].labels(f)
                    << endl;
            }


            if (setPtr)
            {
                setPtr->insert(faceI);
            }

            nErrorPyrs++;
        }

        if (mesh.isInternalFace(faceI))
        {
            // Create the neighbour pyramid - it will have positive volume
            scalar pyrVol =
                pyramidPointFaceRef(f[faceI], cellCentres[nei[faceI]]).mag(p);

            if (pyrVol < minPyrVol)
            {
                if (report)
                {
                    Pout<< "bool polyMeshGeometry::checkFacePyramids("
                        << "const bool, const scalar, const pointField&"
                        << ", const labelList&, labelHashSet*): "
                        << "face " << faceI << " points the wrong way. " << endl
                        << "Pyramid volume: " << -pyrVol
                        << " Face " << f[faceI] << " area: " << f[faceI].mag(p)
                        << " Neighbour cell: " << nei[faceI] << endl
                        << "Neighbour cell vertex labels: "
                        << mesh.cells()[nei[faceI]].labels(f)
                        << endl;
                }

                if (setPtr)
                {
                    setPtr->insert(faceI);
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
            SeriousErrorIn
            (
                "polyMeshGeometry::checkFacePyramids("
                "const bool, const scalar, const pointField&"
                ", const labelList&, labelHashSet*)"
            )   << "Error in face pyramids: faces pointing the wrong way."
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

    for (label faceI = mesh.nInternalFaces(); faceI < mesh.nFaces(); faceI++)
    {
        neiCc[faceI - mesh.nInternalFaces()] = cellCentres[own[faceI]];
    }

    syncTools::swapBoundaryFacePositions(mesh, neiCc);

    label nErrorTets = 0;

    forAll(checkFaces, i)
    {
        label faceI = checkFaces[i];

        // Create the owner pyramid - note: exchange cell and face centre
        // to get positive volume.
        bool tetError = checkFaceTet
        (
            mesh,
            report,
            minTetQuality,
            p,
            faceI,
            cellCentres[own[faceI]],    // face centre
            faceCentres[faceI],         // cell centre
            setPtr
        );

        if (tetError)
        {
            nErrorTets++;
        }

        if (mesh.isInternalFace(faceI))
        {
            // Create the neighbour tets - they will have positive volume
            bool tetError = checkFaceTet
            (
                mesh,
                report,
                minTetQuality,
                p,
                faceI,
                faceCentres[faceI],         // face centre
                cellCentres[nei[faceI]],    // cell centre
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
                    faceI,
                    minTetQuality,
                    report
                ) == -1
            )
            {
                if (setPtr)
                {
                    setPtr->insert(faceI);
                }

                nErrorTets++;
            }
        }
        else
        {
            label patchI = patches.whichPatch(faceI);

            if (patches[patchI].coupled())
            {
                if
                (
                    polyMeshTetDecomposition::findSharedBasePoint
                    (
                        mesh,
                        faceI,
                        neiCc[faceI - mesh.nInternalFaces()],
                        minTetQuality,
                        report
                    ) == -1
                )
                {
                    if (setPtr)
                    {
                        setPtr->insert(faceI);
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
                        faceI,
                        minTetQuality,
                        report
                    ) == -1
                )
                {
                    if (setPtr)
                    {
                        setPtr->insert(faceI);
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
            SeriousErrorIn
            (
                "polyMeshGeometry::checkFaceTets("
                "const bool, const scalar, const pointField&, const pointField&"
                ", const labelList&, labelHashSet*)"
            )   << "Error in face decomposition: negative tets."
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
        label faceI = checkFaces[i];

        if (mesh.isInternalFace(faceI))
        {
            scalar skewness = primitiveMeshTools::faceSkewness
            (
                mesh,
                points,
                faceCentres,
                faceAreas,

                faceI,
                cellCentres[own[faceI]],
                cellCentres[nei[faceI]]
            );

            // Check if the skewness vector is greater than the PN vector.
            // This does not cause trouble but is a good indication of a poor
            // mesh.
            if (skewness > internalSkew)
            {
                if (report)
                {
                    Pout<< "Severe skewness for face " << faceI
                        << " skewness = " << skewness << endl;
                }

                if (setPtr)
                {
                    setPtr->insert(faceI);
                }

                nWarnSkew++;
            }

            maxSkew = max(maxSkew, skewness);
        }
        else if (patches[patches.whichPatch(faceI)].coupled())
        {
            scalar skewness = primitiveMeshTools::faceSkewness
            (
                mesh,
                points,
                faceCentres,
                faceAreas,

                faceI,
                cellCentres[own[faceI]],
                neiCc[faceI-mesh.nInternalFaces()]
            );

            // Check if the skewness vector is greater than the PN vector.
            // This does not cause trouble but is a good indication of a poor
            // mesh.
            if (skewness > internalSkew)
            {
                if (report)
                {
                    Pout<< "Severe skewness for coupled face " << faceI
                        << " skewness = " << skewness << endl;
                }

                if (setPtr)
                {
                    setPtr->insert(faceI);
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

                faceI,
                cellCentres[own[faceI]]
            );


            // Check if the skewness vector is greater than the PN vector.
            // This does not cause trouble but is a good indication of a poor
            // mesh.
            if (skewness > boundarySkew)
            {
                if (report)
                {
                    Pout<< "Severe skewness for boundary face " << faceI
                        << " skewness = " << skewness << endl;
                }

                if (setPtr)
                {
                    setPtr->insert(faceI);
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
            WarningIn
            (
                "polyMeshGeometry::checkFaceSkewness"
                "(const bool, const scalar, const labelList&, labelHashSet*)"
            )   << "Large face skewness detected.  Max skewness = "
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

    for (label faceI = mesh.nInternalFaces(); faceI < mesh.nFaces(); faceI++)
    {
        neiCc[faceI-mesh.nInternalFaces()] = cellCentres[own[faceI]];
    }
    syncTools::swapBoundaryFacePositions(mesh, neiCc);


    scalar minWeight = GREAT;

    label nWarnWeight = 0;

    forAll(checkFaces, i)
    {
        label faceI = checkFaces[i];

        const point& fc = faceCentres[faceI];
        const vector& fa = faceAreas[faceI];

        scalar dOwn = mag(fa & (fc-cellCentres[own[faceI]]));

        if (mesh.isInternalFace(faceI))
        {
            scalar dNei = mag(fa & (cellCentres[nei[faceI]]-fc));
            scalar weight = min(dNei,dOwn)/(dNei+dOwn+VSMALL);

            if (weight < warnWeight)
            {
                if (report)
                {
                    Pout<< "Small weighting factor for face " << faceI
                        << " weight = " << weight << endl;
                }

                if (setPtr)
                {
                    setPtr->insert(faceI);
                }

                nWarnWeight++;
            }

            minWeight = min(minWeight, weight);
        }
        else
        {
            label patchI = patches.whichPatch(faceI);

            if (patches[patchI].coupled())
            {
                scalar dNei = mag(fa & (neiCc[faceI-mesh.nInternalFaces()]-fc));
                scalar weight = min(dNei,dOwn)/(dNei+dOwn+VSMALL);

                if (weight < warnWeight)
                {
                    if (report)
                    {
                        Pout<< "Small weighting factor for face " << faceI
                            << " weight = " << weight << endl;
                    }

                    if (setPtr)
                    {
                        setPtr->insert(faceI);
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
        scalar weight = min(dNei,dOwn)/(dNei+dOwn+VSMALL);

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
            WarningIn
            (
                "polyMeshGeometry::checkFaceWeights"
                "(const bool, const scalar, const labelList&, labelHashSet*)"
            )   << "Small interpolation weight detected.  Min weight = "
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

    for (label faceI = mesh.nInternalFaces(); faceI < mesh.nFaces(); faceI++)
    {
        neiVols[faceI-mesh.nInternalFaces()] = cellVolumes[own[faceI]];
    }
    syncTools::swapBoundaryFaceList(mesh, neiVols);


    scalar minRatio = GREAT;

    label nWarnRatio = 0;

    forAll(checkFaces, i)
    {
        label faceI = checkFaces[i];

        scalar ownVol = mag(cellVolumes[own[faceI]]);

        scalar neiVol = -GREAT;

        if (mesh.isInternalFace(faceI))
        {
            neiVol = mag(cellVolumes[nei[faceI]]);
        }
        else
        {
            label patchI = patches.whichPatch(faceI);

            if (patches[patchI].coupled())
            {
                neiVol = mag(neiVols[faceI-mesh.nInternalFaces()]);
            }
        }

        if (neiVol >= 0)
        {
            scalar ratio = min(ownVol, neiVol) / (max(ownVol, neiVol) + VSMALL);

            if (ratio < warnRatio)
            {
                if (report)
                {
                    Pout<< "Small ratio for face " << faceI
                        << " ratio = " << ratio << endl;
                }

                if (setPtr)
                {
                    setPtr->insert(faceI);
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
            scalar ratio = min(ownVol, neiVol) / (max(ownVol, neiVol) + VSMALL);

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
            WarningIn
            (
                "polyMeshGeometry::checkVolRatio"
                "(const bool, const scalar, const labelList&, labelHashSet*)"
            )   << "Small volume ratio detected.  Min ratio = "
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
    if (maxDeg < -SMALL || maxDeg > 180+SMALL)
    {
        FatalErrorIn
        (
            "polyMeshGeometry::checkFaceAngles"
            "(const bool, const scalar, const pointField&, const labelList&"
            ", labelHashSet*)"
        )   << "maxDeg should be [0..180] but is now " << maxDeg
            << abort(FatalError);
    }

    const scalar maxSin = Foam::sin(degToRad(maxDeg));

    const faceList& fcs = mesh.faces();

    scalar maxEdgeSin = 0.0;

    label nConcave = 0;

    label errorFaceI = -1;

    forAll(checkFaces, i)
    {
        label faceI = checkFaces[i];

        const face& f = fcs[faceI];

        vector faceNormal = faceAreas[faceI];
        faceNormal /= mag(faceNormal) + VSMALL;

        // Get edge from f[0] to f[size-1];
        vector ePrev(p[f.first()] - p[f.last()]);
        scalar magEPrev = mag(ePrev);
        ePrev /= magEPrev + VSMALL;

        forAll(f, fp0)
        {
            // Get vertex after fp
            label fp1 = f.fcIndex(fp0);

            // Normalized vector between two consecutive points
            vector e10(p[f[fp1]] - p[f[fp0]]);
            scalar magE10 = mag(e10);
            e10 /= magE10 + VSMALL;

            if (magEPrev > SMALL && magE10 > SMALL)
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

                    if ((edgeNormal & faceNormal) < SMALL)
                    {
                        if (faceI != errorFaceI)
                        {
                            // Count only one error per face.
                            errorFaceI = faceI;
                            nConcave++;
                        }

                        if (setPtr)
                        {
                            setPtr->insert(faceI);
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
        if (maxEdgeSin > SMALL)
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
            WarningIn
            (
                "polyMeshGeometry::checkFaceAngles"
                "(const bool, const scalar,  const pointField&"
                ", const labelList&, labelHashSet*)"
            )   << nConcave  << " face points with severe concave angle (> "
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
    if (minTwist < -1-SMALL || minTwist > 1+SMALL)
    {
        FatalErrorIn
        (
            "polyMeshGeometry::checkFaceTwist"
            "(const bool, const scalar, const polyMesh&, const pointField&"
            ", const pointField&, const pointField&, const pointField&"
            ", const labelList&, labelHashSet*)"
        )   << "minTwist should be [-1..1] but is now " << minTwist
            << abort(FatalError);
    }


    const faceList& fcs = mesh.faces();


    label nWarped = 0;

//    forAll(checkFaces, i)
//    {
//        label faceI = checkFaces[i];
//
//        const face& f = fcs[faceI];
//
//        scalar magArea = mag(faceAreas[faceI]);
//
//        if (f.size() > 3 && magArea > VSMALL)
//        {
//            const vector nf = faceAreas[faceI] / magArea;
//
//            const point& fc = faceCentres[faceI];
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
//                    ).normal()
//                );
//
//                scalar magTri = mag(triArea);
//
//                if (magTri > VSMALL && ((nf & triArea/magTri) < minTwist))
//                {
//                    nWarped++;
//
//                    if (setPtr)
//                    {
//                        setPtr->insert(faceI);
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

    for (label faceI = mesh.nInternalFaces(); faceI < mesh.nFaces(); faceI++)
    {
        neiCc[faceI-mesh.nInternalFaces()] = cellCentres[own[faceI]];
    }
    syncTools::swapBoundaryFacePositions(mesh, neiCc);

    forAll(checkFaces, i)
    {
        label faceI = checkFaces[i];

        const face& f = fcs[faceI];

        if (f.size() > 3)
        {
            vector nf(vector::zero);

            if (mesh.isInternalFace(faceI))
            {
                nf = cellCentres[nei[faceI]] - cellCentres[own[faceI]];
                nf /= mag(nf) + VSMALL;
            }
            else if (patches[patches.whichPatch(faceI)].coupled())
            {
                nf =
                    neiCc[faceI-mesh.nInternalFaces()]
                  - cellCentres[own[faceI]];
                nf /= mag(nf) + VSMALL;
            }
            else
            {
                nf = faceCentres[faceI] - cellCentres[own[faceI]];
                nf /= mag(nf) + VSMALL;
            }

            if (nf != vector::zero)
            {
                const point& fc = faceCentres[faceI];

                forAll(f, fpI)
                {
                    vector triArea
                    (
                        triPointRef
                        (
                            p[f[fpI]],
                            p[f.nextLabel(fpI)],
                            fc
                        ).normal()
                    );

                    scalar magTri = mag(triArea);

                    if (magTri > VSMALL && ((nf & triArea/magTri) < minTwist))
                    {
                        nWarped++;

                        if (setPtr)
                        {
                            setPtr->insert(faceI);
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
            WarningIn
            (
                "polyMeshGeometry::checkFaceTwist"
                "(const bool, const scalar, const polyMesh&, const pointField&"
                ", const pointField&, const pointField&, const pointField&"
                ", const labelList&, labelHashSet*)"
            )   << nWarped  << " faces with severe warpage "
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
    if (minTwist < -1-SMALL || minTwist > 1+SMALL)
    {
        FatalErrorIn
        (
            "polyMeshGeometry::checkTriangleTwist"
            "(const bool, const scalar, const polyMesh&, const pointField&"
            ", const labelList&, labelHashSet*)"
        )   << "minTwist should be [-1..1] but is now " << minTwist
            << abort(FatalError);
    }

    const faceList& fcs = mesh.faces();

    label nWarped = 0;

    forAll(checkFaces, i)
    {
        label faceI = checkFaces[i];

        const face& f = fcs[faceI];

        if (f.size() > 3)
        {
            const point& fc = faceCentres[faceI];

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
                ).normal();

                scalar magTri = mag(prevN);

                if (magTri > VSMALL)
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
                        ).normal()
                    );
                    scalar magTri = mag(triN);

                    if (magTri > VSMALL)
                    {
                        triN /= magTri;

                        if ((prevN & triN) < minTwist)
                        {
                            nWarped++;

                            if (setPtr)
                            {
                                setPtr->insert(faceI);
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
                            setPtr->insert(faceI);
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
            WarningIn
            (
                "polyMeshGeometry::checkTriangleTwist"
                "(const bool, const scalar, const polyMesh&"
                ", const pointField&, const labelList&, labelHashSet*)"
            )   << nWarped  << " faces with severe warpage "
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
    if (minFlatness < -SMALL || minFlatness > 1+SMALL)
    {
        FatalErrorIn
        (
            "polyMeshGeometry::checkFaceFlatness"
            "(const bool, const scalar, const polyMesh&, const pointField&"
            ", const pointField&, const labelList&, labelHashSet*)"
        )   << "minFlatness should be [0..1] but is now " << minFlatness
            << abort(FatalError);
    }

    const faceList& fcs = mesh.faces();

    label nWarped = 0;

    forAll(checkFaces, i)
    {
        label faceI = checkFaces[i];

        const face& f = fcs[faceI];

        if (f.size() > 3)
        {
            const point& fc = faceCentres[faceI];

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

            if (sumArea/mag(faceAreas[faceI]) < minFlatness)
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

    if (report)
    {
        if (nWarped> 0)
        {
            Info<< "There are " << nWarped
                << " faces with area of invidual triangles"
                << " compared to overall area less than "
                << minFlatness << nl << endl;
        }
        else
        {
            Info<< "All faces are flat in that the area of invidual triangles"
                << " compared to overall area is less than "
                << minFlatness << nl << endl;
        }
    }

    if (nWarped > 0)
    {
        if (report)
        {
            WarningIn
            (
                "polyMeshGeometry::checkFaceFlatness"
                "(const bool, const scalar, const polyMesh&"
                ", const pointField&, const pointField&, const labelList&"
                ", labelHashSet*)"
            )   << nWarped  << " non-flat faces "
                << "(area of invidual triangles"
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
        label faceI = checkFaces[i];

        if (mag(faceAreas[faceI]) < minArea)
        {
            if (setPtr)
            {
                setPtr->insert(faceI);
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
            WarningIn
            (
                "polyMeshGeometry::checkFaceArea"
                "(const bool, const scalar, const polyMesh&"
                ", const pointField&, const labelList&, labelHashSet*)"
            )   << nZeroArea  << " faces with area < " << minArea
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

    scalar minDet = GREAT;
    scalar sumDet = 0.0;
    label nSumDet = 0;
    label nWarnDet = 0;

    forAll(affectedCells, i)
    {
        const cell& cFaces = cells[affectedCells[i]];

        tensor areaSum(tensor::zero);
        scalar magAreaSum = 0;

        forAll(cFaces, cFaceI)
        {
            label faceI = cFaces[cFaceI];

            scalar magArea = mag(faceAreas[faceI]);

            magAreaSum += magArea;
            areaSum += faceAreas[faceI]*(faceAreas[faceI]/(magArea+VSMALL));
        }

        scalar scaledDet = det(areaSum/(magAreaSum+VSMALL))/0.037037037037037;

        minDet = min(minDet, scaledDet);
        sumDet += scaledDet;
        nSumDet++;

        if (scaledDet < warnDet)
        {
            if (setPtr)
            {
                // Insert all faces of the cell.
                forAll(cFaces, cFaceI)
                {
                    label faceI = cFaces[cFaceI];
                    setPtr->insert(faceI);
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
            WarningIn
            (
                "polyMeshGeometry::checkCellDeterminant"
                "(const bool, const scalar, const polyMesh&"
                ", const pointField&, const labelList&, const labelList&"
                ", labelHashSet*)"
            )   << nWarnDet << " cells with determinant < " << warnDet
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
