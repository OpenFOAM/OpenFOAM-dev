/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2014 OpenFOAM Foundation
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
    Foam::processorCyclicPolyPatch

Description
    Neighbour processor patch.

    Note: morph patch face ordering is geometric.

SourceFiles
    processorCyclicPolyPatch.C

\*---------------------------------------------------------------------------*/

#ifndef processorCyclicPolyPatch_H
#define processorCyclicPolyPatch_H

#include "processorPolyPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class processorCyclicPolyPatch Declaration
\*---------------------------------------------------------------------------*/

class processorCyclicPolyPatch
:
    public processorPolyPatch
{
    // Private data

        //- Name of originating patch
        const word referPatchName_;

        //- Message tag to use for communication
        mutable int tag_;

        //- Index of originating patch
        mutable label referPatchID_;


protected:

    // Protected Member functions

        //- Initialise the calculation of the patch geometry
        void initGeometry(PstreamBuffers&);

        //- Calculate the patch geometry
        void calcGeometry(PstreamBuffers&);

        //- Calculate the patch geometry with externally
        //  provided geometry
        virtual void calcGeometry
        (
            const primitivePatch& referPatch,
            const pointField& thisCtrs,
            const vectorField& thisAreas,
            const pointField& thisCc,
            const pointField& nbrCtrs,
            const vectorField& nbrAreas,
            const pointField& nbrCc
        )
        {
            notImplemented("processorCyclicPolyPatch::calcGeometry(..)");
        }

        //- Initialise the patches for moving points
        void initMovePoints(PstreamBuffers&, const pointField&);

        //- Correct patches after moving points
        void movePoints(PstreamBuffers&, const pointField&);

        //- Initialise the update of the patch topology
        virtual void initUpdateMesh(PstreamBuffers&);

        //- Update of the patch topology
        virtual void updateMesh(PstreamBuffers&);


public:

    //- Runtime type information
    TypeName("processorCyclic");


    // Constructors

        //- Construct from components
        processorCyclicPolyPatch
        (
            const word& name,
            const label size,
            const label start,
            const label index,
            const polyBoundaryMesh& bm,
            const int myProcNo,
            const int neighbProcNo,
            const word& referPatchName,
            const transformType transform = UNKNOWN,
            const word& patchType = typeName
        );

        //- Construct from dictionary
        processorCyclicPolyPatch
        (
            const word& name,
            const dictionary& dict,
            const label index,
            const polyBoundaryMesh&,
            const word& patchType
        );

        //- Construct as copy, resetting the boundary mesh
        processorCyclicPolyPatch
        (
            const processorCyclicPolyPatch&,
            const polyBoundaryMesh&
        );

        //- Construct as given the original patch and resetting the
        //  face list and boundary mesh information
        processorCyclicPolyPatch
        (
            const processorCyclicPolyPatch& pp,
            const polyBoundaryMesh& bm,
            const label index,
            const label newSize,
            const label newStart
        );

        //- Construct as given the original patch and resetting the
        //  face list, boundary mesh information and referPatch
        processorCyclicPolyPatch
        (
            const processorCyclicPolyPatch& pp,
            const polyBoundaryMesh& bm,
            const label index,
            const label newSize,
            const label newStart,
            const word& referPatchName
        );

        //- Construct given the original patch and a map
        processorCyclicPolyPatch
        (
            const processorCyclicPolyPatch& pp,
            const polyBoundaryMesh& bm,
            const label index,
            const labelUList& mapAddressing,
            const label newStart
        );


        //- Construct and return a clone, resetting the boundary mesh
        virtual autoPtr<polyPatch> clone(const polyBoundaryMesh& bm) const
        {
            return autoPtr<polyPatch>(new processorCyclicPolyPatch(*this, bm));
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
                new processorCyclicPolyPatch
                (
                    *this,
                    bm,
                    index,
                    newSize,
                    newStart
                )
            );
        }

        //- Construct and return a clone, resetting the face list
        //  and boundary mesh
        virtual autoPtr<polyPatch> clone
        (
            const polyBoundaryMesh& bm,
            const label index,
            const label newSize,
            const label newStart,
            const word& referPatchName
        ) const
        {
            return autoPtr<polyPatch>
            (
                new processorCyclicPolyPatch
                (
                    *this,
                    bm,
                    index,
                    newSize,
                    newStart,
                    referPatchName
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
                new processorCyclicPolyPatch
                (
                    *this,
                    bm,
                    index,
                    mapAddressing,
                    newStart
                )
            );
        }


    // Destructor
    virtual ~processorCyclicPolyPatch();


    // Member functions

        const word& referPatchName() const
        {
            return referPatchName_;
        }

        //- Referring patchID.
        label referPatchID() const
        {
            if (referPatchID_ == -1)
            {
                referPatchID_ = this->boundaryMesh().findPatchID
                (
                    referPatchName_
                );
                if (referPatchID_ == -1)
                {
                    FatalErrorIn
                    (
                        "processorCyclicPolyPatch::referPatchID() const"
                    )   << "Illegal referPatch name " << referPatchName_
                        << endl << "Valid patch names are "
                        << this->boundaryMesh().names()
                        << exit(FatalError);
                }
            }
            return referPatchID_;
        }

        const coupledPolyPatch& referPatch() const
        {
            const polyPatch& pp = this->boundaryMesh()[referPatchID()];
            return refCast<const coupledPolyPatch>(pp);
        }

        //- Return message tag to use for communication
        virtual int tag() const;

        //- Does this side own the patch ?
        virtual bool owner() const
        {
            return referPatch().owner();
        }

        //- Type of transform
        virtual transformType transform() const
        {
            return referPatch().transform();
        }

        //- Type of transform
        //  This is currently only for use when collapsing generated
        //  meshes that can have zero area faces.
        virtual transformType& transform()
        {
            return const_cast<coupledPolyPatch&>(referPatch()).transform();
        }

        //- Transform a patch-based position from other side to this side
        virtual void transformPosition(pointField& l) const
        {
            referPatch().transformPosition(l);
        }

        //- Transform a patch-based position from other side to this side
        virtual void transformPosition(point& l, const label facei) const
        {
            referPatch().transformPosition(l, facei);
        }

        //- Are the planes separated.
        virtual bool separated() const
        {
            return referPatch().separated();
        }

        //- If the planes are separated the separation vector.
        virtual const vectorField& separation() const
        {
            return referPatch().separation();
        }

        //- Are the cyclic planes parallel.
        virtual bool parallel() const
        {
            return referPatch().parallel();
        }

        //- Return face transformation tensor.
        virtual const tensorField& forwardT() const
        {
            return referPatch().forwardT();
        }

        //- Return neighbour-cell transformation tensor.
        virtual const tensorField& reverseT() const
        {
            return referPatch().reverseT();
        }

        //- Are faces collocated. Either size 0,1 or length of patch
        virtual const boolList& collocated() const
        {
            return referPatch().collocated();
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
