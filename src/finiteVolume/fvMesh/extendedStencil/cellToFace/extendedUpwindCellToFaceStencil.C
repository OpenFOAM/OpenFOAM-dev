/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
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

#include "extendedUpwindCellToFaceStencil.H"
#include "cellToFaceStencil.H"
#include "syncTools.H"
#include "SortableList.H"
#include "dummyTransform.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::extendedUpwindCellToFaceStencil::selectOppositeFaces
(
    const boolList& nonEmptyFace,
    const scalar minOpposedness,
    const label facei,
    const label celli,
    DynamicList<label>& oppositeFaces
) const
{
    const vectorField& areas = mesh_.faceAreas();
    const labelList& own = mesh_.faceOwner();
    const cell& cFaces = mesh_.cells()[celli];

    SortableList<scalar> opposedness(cFaces.size(), -great);

    // Pick up all the faces that oppose this one.
    forAll(cFaces, i)
    {
        label otherFacei = cFaces[i];

        if (otherFacei != facei && nonEmptyFace[otherFacei])
        {
            if ((own[otherFacei] == celli) == (own[facei] == celli))
            {
                opposedness[i] = -(areas[otherFacei] & areas[facei]);
            }
            else
            {
                opposedness[i] = (areas[otherFacei] & areas[facei]);
            }
        }
    }

    label sz = opposedness.size();

    oppositeFaces.clear();

    scalar myAreaSqr = magSqr(areas[facei]);

    if (myAreaSqr > vSmall)
    {
        forAll(opposedness, i)
        {
            opposedness[i] /= myAreaSqr;
        }
        // Sort in incrementing order
        opposedness.sort();

        // Pick largest no matter what
        oppositeFaces.append(cFaces[opposedness.indices()[sz-1]]);

        for (label i = sz-2; i >= 0; --i)
        {
            if (opposedness[i] < minOpposedness)
            {
                break;
            }
            oppositeFaces.append(cFaces[opposedness.indices()[i]]);
        }
    }
    else
    {
        // Sort in incrementing order
        opposedness.sort();

        // Tiny face. Do what?
        // Pick largest no matter what
        oppositeFaces.append(cFaces[opposedness.indices()[sz-1]]);
    }
}


void Foam::extendedUpwindCellToFaceStencil::transportStencil
(
    const boolList& nonEmptyFace,
    const labelListList& faceStencil,
    const scalar minOpposedness,
    const label facei,
    const label celli,
    const bool stencilHasNeighbour,

    DynamicList<label>& oppositeFaces,
    labelHashSet& faceStencilSet,
    labelList& transportedStencil
) const
{
    label globalOwn = faceStencil[facei][0];
    label globalNei = -1;
    if (stencilHasNeighbour && faceStencil[facei].size() >= 2)
    {
        globalNei = faceStencil[facei][1];
    }


    selectOppositeFaces
    (
        nonEmptyFace,
        minOpposedness,
        facei,
        celli,
        oppositeFaces
    );

    // Collect all stencils of oppositefaces
    faceStencilSet.clear();
    forAll(oppositeFaces, i)
    {
        const labelList& fStencil = faceStencil[oppositeFaces[i]];

        forAll(fStencil, j)
        {
            label globalI = fStencil[j];

            if (globalI != globalOwn && globalI != globalNei)
            {
                faceStencilSet.insert(globalI);
            }
        }
    }

    // Add my owner and neighbour first.
    if (stencilHasNeighbour)
    {
        transportedStencil.setSize(faceStencilSet.size()+2);
        label n = 0;
        transportedStencil[n++] = globalOwn;
        transportedStencil[n++] = globalNei;

        forAllConstIter(labelHashSet, faceStencilSet, iter)
        {
            if (iter.key() != globalOwn && iter.key() != globalNei)
            {
                transportedStencil[n++] = iter.key();
            }
        }
        if (n != transportedStencil.size())
        {
            FatalErrorInFunction
                << "problem:" << faceStencilSet
                << abort(FatalError);
        }
    }
    else
    {
        transportedStencil.setSize(faceStencilSet.size()+1);
        label n = 0;
        transportedStencil[n++] = globalOwn;

        forAllConstIter(labelHashSet, faceStencilSet, iter)
        {
            if (iter.key() != globalOwn)
            {
                transportedStencil[n++] = iter.key();
            }
        }
        if (n != transportedStencil.size())
        {
            FatalErrorInFunction
                << "problem:" << faceStencilSet
                << abort(FatalError);
        }
    }
}


void Foam::extendedUpwindCellToFaceStencil::transportStencils
(
    const labelListList& faceStencil,
    const scalar minOpposedness,
    labelListList& ownStencil,
    labelListList& neiStencil
)
{
    const polyBoundaryMesh& patches = mesh_.boundaryMesh();
    const label nBnd = mesh_.nFaces()-mesh_.nInternalFaces();
    const labelList& own = mesh_.faceOwner();
    const labelList& nei = mesh_.faceNeighbour();

    // Work arrays
    DynamicList<label> oppositeFaces;
    labelHashSet faceStencilSet;


    // For quick detection of empty faces
    boolList nonEmptyFace(mesh_.nFaces(), true);
    forAll(patches, patchi)
    {
        const polyPatch& pp = patches[patchi];

        if (isA<emptyPolyPatch>(pp))
        {
            label facei = pp.start();
            forAll(pp, i)
            {
                nonEmptyFace[facei++] = false;
            }
        }
    }


    // Do the owner side
    // ~~~~~~~~~~~~~~~~~
    // stencil is synchronised at entry so no need to swap.

    ownStencil.setSize(mesh_.nFaces());

    // Internal faces
    for (label facei = 0; facei < mesh_.nInternalFaces(); facei++)
    {
        // Get stencil as owner + neighbour + stencil from 'opposite' faces
        transportStencil
        (
            nonEmptyFace,
            faceStencil,
            minOpposedness,
            facei,
            own[facei],
            true,                   // stencilHasNeighbour
            oppositeFaces,
            faceStencilSet,
            ownStencil[facei]
        );
    }
    // Boundary faces
    forAll(patches, patchi)
    {
        const polyPatch& pp = patches[patchi];
        label facei = pp.start();

        if (pp.coupled())
        {
            forAll(pp, i)
            {
                transportStencil
                (
                    nonEmptyFace,
                    faceStencil,
                    minOpposedness,
                    facei,
                    own[facei],
                    true,                   // stencilHasNeighbour

                    oppositeFaces,
                    faceStencilSet,
                    ownStencil[facei]
                );
                facei++;
            }
        }
        else if (!isA<emptyPolyPatch>(pp))
        {
            forAll(pp, i)
            {
                // faceStencil does not contain neighbour
                transportStencil
                (
                    nonEmptyFace,
                    faceStencil,
                    minOpposedness,
                    facei,
                    own[facei],
                    false,                  // stencilHasNeighbour

                    oppositeFaces,
                    faceStencilSet,
                    ownStencil[facei]
                );
                facei++;
            }
        }
    }


    // Swap coupled boundary stencil
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    labelListList neiBndStencil(nBnd);
    for (label facei = mesh_.nInternalFaces(); facei < mesh_.nFaces(); facei++)
    {
        neiBndStencil[facei-mesh_.nInternalFaces()] = ownStencil[facei];
    }
    // syncTools::swapBoundaryFaceList(mesh_, neiBndStencil);
    syncTools::syncBoundaryFaceList
    (
        mesh_,
        neiBndStencil,
        eqOp<labelList>(),
        dummyTransform()
    );



    // Do the neighbour side
    // ~~~~~~~~~~~~~~~~~~~~~
    // - internal faces : get opposite faces on neighbour side
    // - boundary faces : empty
    // - coupled faces  : in neiBndStencil

    neiStencil.setSize(mesh_.nFaces());

    // Internal faces
    for (label facei = 0; facei < mesh_.nInternalFaces(); facei++)
    {
        transportStencil
        (
            nonEmptyFace,
            faceStencil,
            minOpposedness,
            facei,
            nei[facei],
            true,                   // stencilHasNeighbour

            oppositeFaces,
            faceStencilSet,
            neiStencil[facei]
        );
    }

    // Boundary faces
    forAll(patches, patchi)
    {
        const polyPatch& pp = patches[patchi];
        label facei = pp.start();

        if (pp.coupled())
        {
            forAll(pp, i)
            {
                neiStencil[facei].transfer
                (
                    neiBndStencil[facei-mesh_.nInternalFaces()]
                );
                facei++;
            }
        }
        else
        {
            // Boundary has empty neighbour stencil
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::extendedUpwindCellToFaceStencil::extendedUpwindCellToFaceStencil
(
    const cellToFaceStencil& stencil,
    const bool pureUpwind,
    const scalar minOpposedness
)
:
    extendedCellToFaceStencil(stencil.mesh()),
    pureUpwind_(pureUpwind)
{
    // forAll(stencil, facei)
    //{
    //    const labelList& fCells = stencil[facei];
    //
    //    Pout<< "Face:" << facei << " at:" << mesh_.faceCentres()[facei]
    //        << endl;
    //
    //    forAll(fCells, i)
    //    {
    //        label globalI = fCells[i];
    //
    //        if (globalI < mesh_.nCells())
    //        {
    //            Pout<< "    cell:" << globalI
    //                << " at:" << mesh_.cellCentres()[globalI] << endl;
    //        }
    //        else
    //        {
    //            label facei = globalI-mesh_.nCells() + mesh_.nInternalFaces();
    //
    //            Pout<< "    boundary:" << facei
    //                << " at:" << mesh_.faceCentres()[facei] << endl;
    //        }
    //    }
    //}
    // Pout<< endl << endl;


    // Transport centred stencil to upwind/downwind face
    transportStencils
    (
        stencil,
        minOpposedness,
        ownStencil_,
        neiStencil_
    );

    {
        List<Map<label>> compactMap(Pstream::nProcs());
        ownMapPtr_.reset
        (
            new mapDistribute
            (
                stencil.globalNumbering(),
                ownStencil_,
                compactMap
            )
        );
    }

    {

        List<Map<label>> compactMap(Pstream::nProcs());
        neiMapPtr_.reset
        (
            new mapDistribute
            (
                stencil.globalNumbering(),
                neiStencil_,
                compactMap
            )
        );
    }

    // stencil now in compact form
    if (pureUpwind_)
    {
        const fvMesh& mesh = dynamic_cast<const fvMesh&>(stencil.mesh());

        List<List<point>> stencilPoints(ownStencil_.size());

        // Owner stencil
        // ~~~~~~~~~~~~~

        collectData(ownMapPtr_(), ownStencil_, mesh.C(), stencilPoints);

        // Mask off all stencil points on wrong side of face
        forAll(stencilPoints, facei)
        {
            const point& fc = mesh.faceCentres()[facei];
            const vector& fArea = mesh.faceAreas()[facei];

            const List<point>& points = stencilPoints[facei];
            const labelList& stencil = ownStencil_[facei];

            DynamicList<label> newStencil(stencil.size());
            forAll(points, i)
            {
                if (((points[i]-fc) & fArea) < 0)
                {
                    newStencil.append(stencil[i]);
                }
            }
            if (newStencil.size() != stencil.size())
            {
                ownStencil_[facei].transfer(newStencil);
            }
        }


        // Neighbour stencil
        // ~~~~~~~~~~~~~~~~~

        collectData(neiMapPtr_(), neiStencil_, mesh.C(), stencilPoints);

        // Mask off all stencil points on wrong side of face
        forAll(stencilPoints, facei)
        {
            const point& fc = mesh.faceCentres()[facei];
            const vector& fArea = mesh.faceAreas()[facei];

            const List<point>& points = stencilPoints[facei];
            const labelList& stencil = neiStencil_[facei];

            DynamicList<label> newStencil(stencil.size());
            forAll(points, i)
            {
                if (((points[i]-fc) & fArea) > 0)
                {
                    newStencil.append(stencil[i]);
                }
            }
            if (newStencil.size() != stencil.size())
            {
                neiStencil_[facei].transfer(newStencil);
            }
        }

        // Note: could compact schedule as well. for if cells are not needed
        // across any boundary anymore. However relatively rare.
    }
}


Foam::extendedUpwindCellToFaceStencil::extendedUpwindCellToFaceStencil
(
    const cellToFaceStencil& stencil
)
:
    extendedCellToFaceStencil(stencil.mesh()),
    pureUpwind_(true)
{
    // Calculate stencil points with full stencil

    ownStencil_ = stencil;

    {
        List<Map<label>> compactMap(Pstream::nProcs());
        ownMapPtr_.reset
        (
            new mapDistribute
            (
                stencil.globalNumbering(),
                ownStencil_,
                compactMap
            )
        );
    }

    const fvMesh& mesh = dynamic_cast<const fvMesh&>(stencil.mesh());

    List<List<point>> stencilPoints(ownStencil_.size());
    collectData(ownMapPtr_(), ownStencil_, mesh.C(), stencilPoints);

    // Split stencil into owner and neighbour
    neiStencil_.setSize(ownStencil_.size());

    forAll(stencilPoints, facei)
    {
        const point& fc = mesh.faceCentres()[facei];
        const vector& fArea = mesh.faceAreas()[facei];

        const List<point>& points = stencilPoints[facei];
        const labelList& stencil = ownStencil_[facei];

        DynamicList<label> newOwnStencil(stencil.size());
        DynamicList<label> newNeiStencil(stencil.size());
        forAll(points, i)
        {
            if (((points[i]-fc) & fArea) > 0)
            {
                newNeiStencil.append(stencil[i]);
            }
            else
            {
                newOwnStencil.append(stencil[i]);
            }
        }
        if (newNeiStencil.size() > 0)
        {
            ownStencil_[facei].transfer(newOwnStencil);
            neiStencil_[facei].transfer(newNeiStencil);
        }
    }

    // Should compact schedule. Or have both return the same schedule.
    neiMapPtr_.reset(new mapDistribute(ownMapPtr_()));
}


// ************************************************************************* //
