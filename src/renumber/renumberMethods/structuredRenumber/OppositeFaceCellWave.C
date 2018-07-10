/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2018 OpenFOAM Foundation
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

#include "OppositeFaceCellWave.H"
#include "polyMesh.H"

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class Type, class TrackingData>
void Foam::OppositeFaceCellWave<Type, TrackingData>::opposingFaceLabels
(
    const label celli,
    const label masterFaceLabel,
    DynamicList<label>& oppositeFaceLabels
) const
{
    // Variant of cell::opposingFaceLabel

    // Algorithm:
    // Go through all the faces of the cell and find the one which
    // does not share a single vertex with the master face.  If there
    // are two or more such faces, return the first one and issue a
    // warning; if there is no opposite face, return -1;

    const face& masterFace = this->mesh_.faces()[masterFaceLabel];

    const labelList& curFaceLabels = this->mesh_.cells()[celli];

    oppositeFaceLabels.clear();

    forAll(curFaceLabels, facei)
    {
        // Compare the face with the master
        const face& curFace = this->mesh_.faces()[curFaceLabels[facei]];

        // Skip the master face
        if (curFaceLabels[facei] != masterFaceLabel)
        {
            bool sharedPoint = false;

            // Compare every vertex of the current face against the
            // vertices of the master face
            forAll(curFace, pointi)
            {
                const label l = curFace[pointi];

                forAll(masterFace, masterPointi)
                {
                    if (masterFace[masterPointi] == l)
                    {
                        sharedPoint = true;
                        break;
                    }
                }

                if (sharedPoint) break;
            }

            // If no points are shared, this is the opposite face
            if (!sharedPoint)
            {
                // Found opposite face
                oppositeFaceLabels.append(curFaceLabels[facei]);
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Iterate, propagating changedFacesInfo across mesh, until no change (or
// maxIter reached). Initial cell values specified.
template<class Type, class TrackingData>
Foam::OppositeFaceCellWave<Type, TrackingData>::OppositeFaceCellWave
(
    const polyMesh& mesh,
    const labelList& changedFaces,
    const List<Type>& changedFacesInfo,
    UList<Type>& allFaceInfo,
    UList<Type>& allCellInfo,
    const label maxIter,
    TrackingData& td
)
:
    FaceCellWave<Type, TrackingData>
    (
        mesh,
        changedFaces,
        changedFacesInfo,
        allFaceInfo,
        allCellInfo,
        0,              // maxIter,
        td
    ),
    changedOppositeFaces_(this->mesh_.nCells())
{
    // Iterate until nothing changes
    label iter = this->iterate(maxIter);

    if ((maxIter > 0) && (iter >= maxIter))
    {
        FatalErrorInFunction
            << "Maximum number of iterations reached. Increase maxIter."
            << endl
            << "    maxIter:" << maxIter << endl
            << "    nChangedCells:" << this->changedCells_.size() << endl
            << "    nChangedFaces:" << this->changedFaces_.size() << endl
            << exit(FatalError);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type, class TrackingData>
Foam::label Foam::OppositeFaceCellWave<Type, TrackingData>::faceToCell()
{
    const labelList& owner = this->mesh_.faceOwner();
    const labelList& neighbour = this->mesh_.faceNeighbour();
    label nInternalFaces = this->mesh_.nInternalFaces();

    DynamicList<label> oppositeFaceLabels;

    forAll(this->changedFaces_, changedFacei)
    {
        label facei = this->changedFaces_[changedFacei];

        if (!this->changedFace_[facei])
        {
            FatalErrorInFunction
                << "Face " << facei
                << " not marked as having been changed"
                << abort(FatalError);
        }


        const Type& neighbourWallInfo = this->allFaceInfo_[facei];

        // Evaluate all connected cells

        // Owner
        {
            label celli = owner[facei];
            Type& currentWallInfo = this->allCellInfo_[celli];

            if (!currentWallInfo.equal(neighbourWallInfo, this->td_))
            {
                // Check if cell is prismatic w.r.t facei
                opposingFaceLabels(celli, facei, oppositeFaceLabels);

                if (oppositeFaceLabels.size())
                {
                    label sz = this->changedCells_.size();
                    this->updateCell
                    (
                        celli,
                        facei,
                        neighbourWallInfo,
                        this->propagationTol_,
                        currentWallInfo
                    );
                    if (this->changedCells_.size() > sz)
                    {
                        label oppFacei = -1;
                        if (oppositeFaceLabels.size() == 1)
                        {
                            oppFacei = oppositeFaceLabels[0];
                        }
                        changedOppositeFaces_.append(oppFacei);
                    }
                }
            }
        }

        // Neighbour.
        if (facei < nInternalFaces)
        {
            label celli = neighbour[facei];
            Type& currentWallInfo2 = this->allCellInfo_[celli];

            if (!currentWallInfo2.equal(neighbourWallInfo, this->td_))
            {
                // Check if cell is prismatic w.r.t facei
                opposingFaceLabels(celli, facei, oppositeFaceLabels);

                if (oppositeFaceLabels.size())
                {
                    label sz = this->changedCells_.size();
                    this->updateCell
                    (
                        celli,
                        facei,
                        neighbourWallInfo,
                        this->propagationTol_,
                        currentWallInfo2
                    );
                    if (this->changedCells_.size() > sz)
                    {
                        label oppFacei = -1;
                        if (oppositeFaceLabels.size() == 1)
                        {
                            oppFacei = oppositeFaceLabels[0];
                        }
                        changedOppositeFaces_.append(oppFacei);
                    }
                }
            }
        }

        // Reset status of face
        this->changedFace_[facei] = false;
    }

    // Handled all changed faces by now
    this->changedFaces_.clear();

    if (debug & 2)
    {
        Pout<< " Changed cells            : " << this->changedCells_.size()
            << endl;
    }

    // Sum changedCells over all procs
    label totNChanged = this->changedCells_.size();

    reduce(totNChanged, sumOp<label>());

    return totNChanged;
}


template<class Type, class TrackingData>
Foam::label Foam::OppositeFaceCellWave<Type, TrackingData>::cellToFace()
{
    forAll(this->changedCells_, changedCelli)
    {
        label celli = this->changedCells_[changedCelli];
        label facei = changedOppositeFaces_[changedCelli];

        if (!this->changedCell_[celli])
        {
            FatalErrorInFunction
                << "Cell " << celli << " not marked as having been changed"
                << abort(FatalError);
        }

        if (facei != -1)
        {
            const Type& neighbourWallInfo = this->allCellInfo_[celli];

            // Evaluate facei

            Type& currentWallInfo = this->allFaceInfo_[facei];

            if (!currentWallInfo.equal(neighbourWallInfo, this->td_))
            {
                this->updateFace
                (
                    facei,
                    celli,
                    neighbourWallInfo,
                    this->propagationTol_,
                    currentWallInfo
                );
            }
        }

        // Reset status of cell
        this->changedCell_[celli] = false;
    }

    // Handled all changed cells by now
    this->changedCells_.clear();
    changedOppositeFaces_.clear();

    if (this->hasCyclicPatches_)
    {
        // Transfer changed faces across cyclic halves
        this->handleCyclicPatches();
    }

    if (this->hasCyclicAMIPatches_)
    {
        this->handleAMICyclicPatches();
    }

    if (Pstream::parRun())
    {
        // Transfer changed faces from neighbouring processors.
        this->handleProcPatches();
    }

    if (debug & 2)
    {
        Pout<< " Changed faces            : " << this->changedFaces_.size()
            << endl;
    }

    // Sum nChangedFaces over all procs
    label totNChanged = this->changedFaces_.size();

    reduce(totNChanged, sumOp<label>());

    return totNChanged;
}


// ************************************************************************* //
