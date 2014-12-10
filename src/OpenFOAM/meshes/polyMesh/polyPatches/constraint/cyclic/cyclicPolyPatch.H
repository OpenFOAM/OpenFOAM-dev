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
    public coupledPolyPatch
{
    // Private data

        //- Name of other half
        mutable word neighbPatchName_;

        //- Optional patchGroup to find neighbPatch
        const coupleGroupIdentifier coupleGroup_;

        //- Index of other half
        mutable label neighbPatchID_;

        // For rotation

            //- Axis of rotation for rotational cyclics
            vector rotationAxis_;

            //- point on axis of rotation for rotational cyclics
            point rotationCentre_;

        // For translation

            //- Translation vector
            vector separationVector_;


        //- List of edges formed from connected points. e[0] is the point on
        //  the first half of the patch, e[1] the corresponding point on the
        //  second half.
        mutable edgeList* coupledPointsPtr_;

        //- List of connected edges. e[0] is the edge on the first half of the
        //  patch, e[1] the corresponding edge on the second half.
        mutable edgeList* coupledEdgesPtr_;

        //- Temporary storage of owner side patch during ordering.
        mutable autoPtr<primitivePatch> ownerPatchPtr_;


    // Private Member Functions

        //- Find amongst selected faces the one with the largest area
        static label findMaxArea(const pointField&, const faceList&);

        void calcTransforms
        (
            const primitivePatch& half0,
            const pointField& half0Ctrs,
            const vectorField& half0Areas,
            const pointField& half1Ctrs,
            const vectorField& half1Areas
        );

        // Face ordering

            // Given a split of faces into left and right half calculate the
            // centres and anchor points. Transform the left points so they
            // align with the right ones
            void getCentresAndAnchors
            (
                const primitivePatch& pp0,
                const primitivePatch& pp1,

                pointField& half0Ctrs,
                pointField& half1Ctrs,
                pointField& anchors0,
                scalarField& tols
            ) const;

            //- Return normal of face at max distance from rotation axis
            vector findFaceMaxRadius(const pointField& faceCentres) const;


protected:

    // Protected Member functions

        //- Recalculate the transformation tensors
        virtual void calcTransforms();

        //- Initialise the calculation of the patch geometry
        virtual void initGeometry(PstreamBuffers&);

        //- Initialise the calculation of the patch geometry
        virtual void initGeometry
        (
            const primitivePatch& referPatch,
            pointField& nbrCtrs,
            vectorField& nbrAreas,
            pointField& nbrCc
        );

        //- Calculate the patch geometry
        virtual void calcGeometry(PstreamBuffers&);

        //- Calculate the patch geometry
        virtual void calcGeometry
        (
            const primitivePatch& referPatch,
            const pointField& thisCtrs,
            const vectorField& thisAreas,
            const pointField& thisCc,
            const pointField& nbrCtrs,
            const vectorField& nbrAreas,
            const pointField& nbrCc
        );

        //- Initialise the patches for moving points
        virtual void initMovePoints(PstreamBuffers&, const pointField&);

        //- Correct patches after moving points
        virtual void movePoints(PstreamBuffers&, const pointField&);

        //- Initialise the update of the patch topology
        virtual void initUpdateMesh(PstreamBuffers&);

        //- Update of the patch topology
        virtual void updateMesh(PstreamBuffers&);

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
            const word& patchType,
            const transformType transform = UNKNOWN
        );

        //- Construct from components
        cyclicPolyPatch
        (
            const word& name,
            const label size,
            const label start,
            const label index,
            const polyBoundaryMesh& bm,
            const word& neighbPatchName,
            const transformType transform,  // transformation type
            const vector& rotationAxis,     // for rotation only
            const point& rotationCentre,    // for rotation only
            const vector& separationVector  // for translation only
        );

        //- Construct from dictionary
        cyclicPolyPatch
        (
            const word& name,
            const dictionary& dict,
            const label index,
            const polyBoundaryMesh& bm,
            const word& patchType
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
            const word& neighbPatchName
        );

        //- Construct given the original patch and a map
        cyclicPolyPatch
        (
            const cyclicPolyPatch& pp,
            const polyBoundaryMesh& bm,
            const label index,
            const labelUList& mapAddressing,
            const label newStart
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
                    neighbPatchName_
                )
            );
        }

        //- Construct and return a clone, resetting the face list
        //  and boundary mesh
        virtual autoPtr<polyPatch> clone
        (
            const polyBoundaryMesh& bm,
            const label index,
            const labelUList& mapAddressing,
            const label newStart
        ) const
        {
            return autoPtr<polyPatch>
            (
                new cyclicPolyPatch(*this, bm, index, mapAddressing, newStart)
            );
        }


    //- Destructor
    virtual ~cyclicPolyPatch();


    // Member Functions

        //- Neighbour patch name
        const word& neighbPatchName() const;

        //- Neighbour patchID
        virtual label neighbPatchID() const;

        virtual bool owner() const
        {
            return index() < neighbPatchID();
        }

        virtual bool neighbour() const
        {
            return !owner();
        }

        const cyclicPolyPatch& neighbPatch() const
        {
            const polyPatch& pp = this->boundaryMesh()[neighbPatchID()];
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

        //- Transform a patch-based position from other side to this side
        virtual void transformPosition(pointField& l) const;

        //- Transform a patch-based position from other side to this side
        virtual void transformPosition(point&, const label facei) const;


        // Transformation

        label transformGlobalFace(const label facei) const
        {
            label offset = facei-start();
            label neighbStart = neighbPatch().start();

            if (offset >= 0 && offset < size())
            {
                return neighbStart+offset;
            }
            else
            {
                FatalErrorIn
                (
                    "cyclicPolyPatch::transformGlobalFace(const label) const"
                )   << "Face " << facei << " not in patch " << name()
                    << exit(FatalError);
                return -1;
            }
        }

        //- Axis of rotation for rotational cyclics
        const vector& rotationAxis() const
        {
            return rotationAxis_;
        }

        //- point on axis of rotation for rotational cyclics
        const point& rotationCentre() const
        {
            return rotationCentre_;
        }

        //- Translation vector for translational cyclics
        const vector& separationVector() const
        {
            return separationVector_;
        }


        //- Initialize ordering for primitivePatch. Does not
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
