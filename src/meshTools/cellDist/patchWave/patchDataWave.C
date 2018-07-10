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

#include "patchDataWave.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// Set initial set of changed faces (= all wall faces)
template<class TransferType>
void Foam::patchDataWave<TransferType>::setChangedFaces
(
    const labelHashSet& patchIDs,
    labelList& changedFaces,
    List<TransferType>& faceDist
) const
{
    const polyMesh& mesh = cellDistFuncs::mesh();

    label nChangedFaces = 0;

    forAll(mesh.boundaryMesh(), patchi)
    {
        if (patchIDs.found(patchi))
        {
            const polyPatch& patch = mesh.boundaryMesh()[patchi];

            const Field<Type>& patchField = initialPatchValuePtrs_[patchi];

            forAll(patch.faceCentres(), patchFacei)
            {
                label meshFacei = patch.start() + patchFacei;

                changedFaces[nChangedFaces] = meshFacei;

                faceDist[nChangedFaces] =
                    TransferType
                    (
                        patch.faceCentres()[patchFacei],
                        patchField[patchFacei],
                        0.0
                    );

                nChangedFaces++;
            }
        }
    }
}


// Copy from MeshWave data into *this (distance) and field_ (transported data)
template<class TransferType>
Foam::label Foam::patchDataWave<TransferType>::getValues
(
    const MeshWave<TransferType>& waveInfo
)
{
    const polyMesh& mesh = cellDistFuncs::mesh();

    const List<TransferType>& cellInfo = waveInfo.allCellInfo();
    const List<TransferType>& faceInfo = waveInfo.allFaceInfo();

    label nIllegal = 0;

    // Copy cell values
    distance_.setSize(cellInfo.size());

    forAll(cellInfo, celli)
    {
        const TransferType & wpn = cellInfo[celli];

        scalar dist = wpn.distSqr();

        if (cellInfo[celli].valid(waveInfo.data()))
        {
            distance_[celli] = Foam::sqrt(dist);

            cellData_[celli] = cellInfo[celli].data();
        }
        else
        {
            // Illegal/unset value. What to do with data?
            // Note: mag for now. Should maybe be member of TransferType?

            distance_[celli] = mag(dist);

            // cellData_[celli] = point::max;
            cellData_[celli] = cellInfo[celli].data();

            nIllegal++;
        }
    }

    // Copy boundary values
    forAll(patchDistance_, patchi)
    {
        const polyPatch& patch = mesh.boundaryMesh()[patchi];

        // Allocate storage for patchDistance
        scalarField* patchFieldPtr = new scalarField(patch.size());

        patchDistance_.set(patchi, patchFieldPtr);

        scalarField& patchField = *patchFieldPtr;

        // Allocate storage for patchData
        Field<Type>* patchDataFieldPtr = new Field<Type>(patch.size());

        patchData_.set(patchi, patchDataFieldPtr);

        Field<Type>& patchDataField = *patchDataFieldPtr;

        // Copy distance and data
        forAll(patchField, patchFacei)
        {
            label meshFacei = patch.start() + patchFacei;

            scalar dist = faceInfo[meshFacei].distSqr();

            if (faceInfo[meshFacei].valid(waveInfo.data()))
            {
                // Adding small to avoid problems with /0 in the turbulence
                // models
                patchField[patchFacei] = Foam::sqrt(dist) + small;

                patchDataField[patchFacei] = faceInfo[meshFacei].data();
            }
            else
            {
                // Illegal/unset value. What to do with data?

                patchField[patchFacei] = mag(dist);

                // patchDataField[patchFacei] = point::max;
                patchDataField[patchFacei] = faceInfo[meshFacei].data();

                nIllegal++;
            }
        }
    }

    return nIllegal;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
template<class TransferType>
Foam::patchDataWave<TransferType>::patchDataWave
(
    const polyMesh& mesh,
    const labelHashSet& patchIDs,
    const UPtrList<Field<Type>>& initialPatchValuePtrs,
    const bool correctWalls
)
:
    cellDistFuncs(mesh),
    patchIDs_(patchIDs),
    initialPatchValuePtrs_(initialPatchValuePtrs),
    correctWalls_(correctWalls),
    nUnset_(0),
    distance_(mesh.nCells()),
    patchDistance_(mesh.boundaryMesh().size()),
    cellData_(mesh.nCells()),
    patchData_(mesh.boundaryMesh().size())
{
    patchDataWave<TransferType>::correct();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class TransferType>
Foam::patchDataWave<TransferType>::~patchDataWave()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Correct for mesh geom/topo changes
template<class TransferType>
void Foam::patchDataWave<TransferType>::correct()
{
    //
    // Set initial changed faces: set TransferType for wall faces
    // to wall centre.
    //

    // Count walls
    label nWalls = sumPatchSize(patchIDs_);

    List<TransferType> faceDist(nWalls);
    labelList changedFaces(nWalls);

    setChangedFaces(patchIDs_, changedFaces, faceDist);

    //
    // Do calculate wall distance by 'growing' from faces.
    //

    MeshWave<TransferType> waveInfo
    (
        mesh(),
        changedFaces,
        faceDist,
        mesh().globalData().nTotalCells()+1 // max iterations
    );


    //
    // Copy distance into return field
    //

    nUnset_ = getValues(waveInfo);

    //
    // Correct wall cells for true distance
    //

    if (correctWalls_)
    {
        Map<label> nearestFace(2 * nWalls);

        // Get distance and indices of nearest face
        correctBoundaryFaceCells
        (
            patchIDs_,
            distance_,
            nearestFace
        );

        correctBoundaryPointCells
        (
            patchIDs_,
            distance_,
            nearestFace
        );

        // Transfer data from nearest face to cell
        const List<TransferType>& faceInfo = waveInfo.allFaceInfo();

        const labelList wallCells(nearestFace.toc());

        forAll(wallCells, wallCelli)
        {
            label celli = wallCells[wallCelli];

            label facei = nearestFace[celli];

            cellData_[celli] = faceInfo[facei].data();
        }
    }
}


// ************************************************************************* //
