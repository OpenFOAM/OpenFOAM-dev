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
    Foam::coupledPolyPatch

Description
    The coupledPolyPatch is an abstract base class for patches that couple
    regions of the computational domain e.g. cyclic and processor-processor
    links.

SourceFiles
    coupledPolyPatch.C

\*---------------------------------------------------------------------------*/

#ifndef coupledPolyPatch_H
#define coupledPolyPatch_H

#include "polyPatch.H"
#include "transformer.H"
#include "diagTensorField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class coupledPolyPatch Declaration
\*---------------------------------------------------------------------------*/

class coupledPolyPatch
:
    public polyPatch
{
private:

    // Private Member Functions

        //- Walk the given primitive patch. Starts at the given seed face and
        //  point within that face. Traverses through the next available
        //  manifold edge to an as yet unvisited face. Stores the order in
        //  which faces are visited and the point/edge through which they are
        //  first walked into. Direction reverses the order in which face-edges
        //  are considered. For two conformal patches, with the same seed face
        //  and point, the maps are guaranteed to match. They can therefore be
        //  used for ordering purposes. This can also be used to identify
        //  contiguous regions of the patch.
        bool walk
        (
            const primitivePatch& pp,
            const bool direction,
            const label seedFacei,
            const label seedFacePointi,
            labelList& faceMap,
            labelList& facePointMap,
            label& mapFacei,
            autoPtr<labelListList>& walks
        ) const;


protected:

    // Protected Classes

        //- Data to pass from owner.initOrder to owner.order
        struct ownToOwnOrderData
        {
            labelList seedFaceis;
        };

        //- Data to pass from owner.initOrder to nbr.order
        struct ownToNbrOrderData
        {
            List<pointField> seedFacePoints;

            inline void transform(const transformer& tr)
            {
                forAll(seedFacePoints, regioni)
                {
                    tr.transformPosition
                    (
                        seedFacePoints[regioni],
                        seedFacePoints[regioni]
                    );
                }
            }
        };

        //- Data to pass from owner.initOrder to nbr.order if debugging
        struct ownToNbrDebugOrderData
        {
            label nFaces;
            label nPoints;
            label nEdges;
            label nInternalEdges;

            inline void transform(const transformer& tr)
            {}
        };


    // Protected Data

        //- Default matching tolerance
        static const scalar defaultMatchTol_;

        //- Local matching tolerance
        const scalar matchTolerance_;

        //- Data to pass from owner.initOrder to owner.order
        mutable autoPtr<ownToOwnOrderData> ownToOwnOrderDataPtr_;


    // Protected Member Functions

        //- Initialise the calculation of the patch geometry
        virtual void initCalcGeometry(PstreamBuffers&) = 0;

        //- Calculate the patch geometry
        virtual void calcGeometry(PstreamBuffers&) = 0;

        //- Initialise the patches for moving points
        virtual void initMovePoints(PstreamBuffers&, const pointField&) = 0;

        //- Correct patches after moving points
        virtual void movePoints(PstreamBuffers&, const pointField&) = 0;

        //- Initialise the update of the patch topology
        virtual void initTopoChange(PstreamBuffers&) = 0;

        //- Update of the patch topology
        virtual void topoChange(PstreamBuffers&) = 0;

        //- Write a patch in OBJ format
        static void writeOBJ
        (
            const fileName&,
            const primitivePatch&
        );

        //- Write lines between two lists of points in OBJ format
        static void writeOBJ
        (
            const fileName&,
            const pointField&,
            const pointField&
        );

        //- Write a set of paths in OBJ format
        static void writeOBJ
        (
            const fileName&,
            const pointField&,
            const labelListList&
        );

        //- Initialise ordering for the given primitivePatch. Fills the
        //  referenced data structures, but leaves transferring them to the
        //  opposite patch to the caller.
        virtual void initOrder
        (
            ownToNbrOrderData& ownToNbr,
            autoPtr<ownToNbrDebugOrderData>& ownToNbrDebugPtr,
            const primitivePatch&
        ) const;

        //- Return new ordering for the given primitivePatch.
        //  Ordering is -faceMap: for every face
        //  index of the new face -rotation:for every new face the clockwise
        //  shift of the original face. Return false if nothing changes
        //  (faceMap is identity, rotation is 0), true otherwise.
        virtual bool order
        (
            const ownToNbrOrderData& ownToNbr,
            const autoPtr<ownToNbrDebugOrderData>& ownToNbrDebugPtr,
            const primitivePatch&,
            labelList& faceMap,
            labelList& rotation
        ) const;


public:

    //- Runtime type information
    TypeName("coupled");


    // Constructors

        //- Construct from components
        coupledPolyPatch
        (
            const word& name,
            const label size,
            const label start,
            const label index,
            const polyBoundaryMesh& bm,
            const word& patchType
        );

        //- Construct from dictionary
        coupledPolyPatch
        (
            const word& name,
            const dictionary& dict,
            const label index,
            const polyBoundaryMesh& bm,
            const word& patchType
        );

        //- Construct as copy, resetting the boundary mesh
        coupledPolyPatch(const coupledPolyPatch&, const polyBoundaryMesh&);

        //- Construct given the original patch and resetting the
        //  face list and boundary mesh information
        coupledPolyPatch
        (
            const coupledPolyPatch& pp,
            const polyBoundaryMesh& bm,
            const label index,
            const label newSize,
            const label newStart
        );


    //- Destructor
    virtual ~coupledPolyPatch();


    // Member Functions

        // Access

            //- Return true because this patch is coupled
            virtual bool coupled() const
            {
                return true;
            }

            //- Does this side own the patch ?
            virtual bool owner() const = 0;

            //- Does the coupled side own the patch ?
            virtual bool neighbour() const
            {
                return !owner();
            }

            //- Return transformation between the coupled patches
            virtual const transformer& transform() const = 0;

            //- Return the matching tolerance
            scalar matchTolerance() const
            {
                return matchTolerance_;
            }


        //- Initialise ordering for primitivePatch. Does not
        //  refer to *this (except for name() and type() etc.)
        virtual void initOrder
        (
            PstreamBuffers&,
            const primitivePatch&
        ) const = 0;

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
        ) const = 0;

        //- Calculate typical tolerance per face. Is currently max distance
        //  from face centre to any of the face vertices.
        static scalarField calcFaceTol
        (
            const UList<face>& faces,
            const pointField& points,
            const pointField& faceCentres
        );

        //- Write the polyPatch data as a dictionary
        virtual void write(Ostream&) const;


    // IOstream Operators

        friend Istream& operator>>(Istream&, ownToNbrOrderData&);
        friend Ostream& operator<<(Ostream&, const ownToNbrOrderData&);
        friend Istream& operator>>(Istream&, ownToNbrDebugOrderData&);
        friend Ostream& operator<<(Ostream&, const ownToNbrDebugOrderData&);
};


Istream& operator>>(Istream&, coupledPolyPatch::ownToNbrOrderData&);
Ostream& operator<<(Ostream&, const coupledPolyPatch::ownToNbrOrderData&);
Istream& operator>>(Istream&, coupledPolyPatch::ownToNbrDebugOrderData&);
Ostream& operator<<(Ostream&, const coupledPolyPatch::ownToNbrDebugOrderData&);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
