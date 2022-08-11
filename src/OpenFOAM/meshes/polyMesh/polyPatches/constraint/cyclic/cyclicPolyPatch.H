/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2022 OpenFOAM Foundation
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

Class
    Foam::cyclicPolyPatch

Description
    Cyclic plane patch.

    Note: morph patch face ordering uses geometric matching so with the
    following restrictions:
        -coupled patches should be flat planes.
        -no rotation in patch plane

    Uses coupledPolyPatch::calcFaceTol to calculate
    tolerance per face which might need tweaking.

    Switch on 'cyclicPolyPatch' debug flag to write .obj files to show
    the matching.

SourceFiles
    cyclicPolyPatch.C

\*---------------------------------------------------------------------------*/

#ifndef cyclicPolyPatch_H
#define cyclicPolyPatch_H

#include "coupledPolyPatch.H"
#include "cyclicTransform.H"
#include "edgeList.H"
#include "polyBoundaryMesh.H"
#include "diagTensorField.H"
#include "coupleGroupIdentifier.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class cyclicPolyPatch Declaration
\*---------------------------------------------------------------------------*/

class cyclicPolyPatch
:
    public coupledPolyPatch,
    public cyclicTransform
{
    // Private Classes

        //- Data to pass from owner.initOrder to nbr.order
        struct ownToNbrCyclicOrderData
        {
            point ctr;
            vector area;
        };


    // Private Data

        //- Name of cyclic neighbour patch
        mutable word nbrPatchName_;

        //- Optional patchGroup to find neighbour patch
        const coupleGroupIdentifier coupleGroup_;

        //- Index of cyclic neighbour patch
        mutable label nbrPatchID_;

        //- List of edges formed from connected points. e[0] is the point on
        //  this patch, e[1] the corresponding point on the cyclic neighbour.
        mutable edgeList* coupledPointsPtr_;

        //- List of connected edges. e[0] is the edge on this patch,
        //  e[1] the corresponding edge on the cyclic neighbour
        mutable edgeList* coupledEdgesPtr_;

        //- Data to pass from owner.initOrder to nbr.order
        mutable autoPtr<ownToNbrOrderData> ownToNbrOrderDataPtr_;

        //- Data to pass from owner.initOrder to nbr.order
        mutable autoPtr<ownToNbrCyclicOrderData> ownToNbrCyclicOrderDataPtr_;

        //- Data to pass from owner.initOrder to nbr.order if debugging
        mutable autoPtr<ownToNbrDebugOrderData> ownToNbrDebugOrderDataPtr_;


protected:

    // Protected Member functions

        //- Initialise the calculation of the patch geometry
        virtual void initCalcGeometry(PstreamBuffers&);

        //- Initialise the calculation of the patch geometry
        virtual void initCalcGeometry
        (
            const primitivePatch& referPatch,
            pointField& nbrCtrs,
            vectorField& nbrAreas,
            pointField& nbrCc
        );

        //- Calculate the patch geometry
        virtual void calcGeometry(PstreamBuffers&);

        //- Initialise the patches for moving points
        virtual void initMovePoints(PstreamBuffers&, const pointField&);

        //- Correct patches after moving points
        virtual void movePoints(PstreamBuffers&, const pointField&);

        //- Initialise the update of the patch topology
        virtual void initTopoChange(PstreamBuffers&);

        //- Update of the patch topology
        virtual void topoChange(PstreamBuffers&);

        //- Reset the patch name
        virtual void rename(const wordList& newNames);

        //- Reset the patch index
        virtual void reorder(const labelUList& newToOldIndex);


public:

    //- Declare friendship with processorCyclicPolyPatch
    friend class processorCyclicPolyPatch;


    //- Runtime type information
    TypeName("cyclic");


    // Constructors

        //- Construct from components
        cyclicPolyPatch
        (
            const word& name,
            const label size,
            const label start,
            const label index,
            const polyBoundaryMesh& bm,
            const word& patchType
        );

        //- Construct from components
        cyclicPolyPatch
        (
            const word& name,
            const label size,
            const label start,
            const label index,
            const polyBoundaryMesh& bm,
            const word& patchType,
            const word& nbrPatchName,
            const cyclicTransform& transform=cyclicTransform(false)
        );

        //- Construct from dictionary
        cyclicPolyPatch
        (
            const word& name,
            const dictionary& dict,
            const label index,
            const polyBoundaryMesh& bm,
            const word& patchType,
            const bool cyclicTransformDefaultIsNone=false
        );

        //- Construct as copy, resetting the boundary mesh
        cyclicPolyPatch(const cyclicPolyPatch&, const polyBoundaryMesh&);

        //- Construct given the original patch and resetting the
        //  face list and boundary mesh information
        cyclicPolyPatch
        (
            const cyclicPolyPatch& pp,
            const polyBoundaryMesh& bm,
            const label index,
            const label newSize,
            const label newStart,
            const word& nbrPatchName
        );

        //- Construct and return a clone, resetting the boundary mesh
        virtual autoPtr<polyPatch> clone(const polyBoundaryMesh& bm) const
        {
            return autoPtr<polyPatch>(new cyclicPolyPatch(*this, bm));
        }

        //- Construct and return a clone, resetting the face list
        //  and boundary mesh
        virtual autoPtr<polyPatch> clone
        (
            const polyBoundaryMesh& bm,
            const label index,
            const label newSize,
            const label newStart
        ) const
        {
            return autoPtr<polyPatch>
            (
                new cyclicPolyPatch
                (
                    *this,
                    bm,
                    index,
                    newSize,
                    newStart,
                    nbrPatchName_
                )
            );
        }


    //- Destructor
    virtual ~cyclicPolyPatch();


    // Member Functions

        //- Neighbour patch name
        const word& nbrPatchName() const;

        //- Neighbour patchID
        virtual label nbrPatchID() const;

        virtual bool owner() const
        {
            return index() < nbrPatchID();
        }

        virtual bool neighbour() const
        {
            return !owner();
        }

        const cyclicPolyPatch& nbrPatch() const
        {
            const polyPatch& pp = this->boundaryMesh()[nbrPatchID()];
            return refCast<const cyclicPolyPatch>(pp);
        }

        //- Return connected points (from patch local to neighbour patch local)
        //  Demand driven calculation. Does primitivePatch::clearOut after
        //  calculation!
        const edgeList& coupledPoints() const;

        //- Return connected edges (from patch local to neighbour patch local).
        //  Demand driven calculation. Does primitivePatch::clearOut after
        //  calculation!
        const edgeList& coupledEdges() const;

        //- Return transformation between the coupled patches
        virtual const transformer& transform() const
        {
            return cyclicTransform::transform();
        }

        //- For a given patch face index, return the corresponding index of the
        //  face on the neighbour
        label transformGlobalFace(const label facei) const
        {
            label offset = facei - start();
            label neiStart = nbrPatch().start();

            if (offset >= 0 && offset < size())
            {
                return neiStart + offset;
            }
            else
            {
                FatalErrorInFunction
                    << "Face " << facei << " not in patch " << name()
                    << exit(FatalError);
                return -1;
            }
        }

        //- Initialise ordering for primitivePatch. Does not
        //  refer to *this (except for name() and type() etc.)
        virtual void initOrder(PstreamBuffers&, const primitivePatch&) const;

        //- Return new ordering for primitivePatch.
        //  Ordering is -faceMap: for every face
        //  index of the new face -rotation:for every new face the clockwise
        //  shift of the original face. Return false if nothing changes
        //  (faceMap is identity, rotation is 0), true otherwise.
        virtual bool order
        (
            PstreamBuffers&,
            const primitivePatch&,
            labelList& faceMap,
            labelList& rotation
        ) const;

        //- Write the polyPatch data as a dictionary
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
