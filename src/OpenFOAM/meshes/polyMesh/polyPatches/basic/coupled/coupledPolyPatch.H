/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
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
public:

    enum transformType
    {
        UNKNOWN,            // unspecified; automatic ordering
        ROTATIONAL,         // rotation along coordinate axis
        TRANSLATIONAL,      // translation
        COINCIDENTFULLMATCH,// assume no transforms
                            // and check the points in faces match
        NOORDERING          // unspecified, no automatic ordering
    };

    static const NamedEnum<transformType, 5> transformTypeNames;


private:

    // Private Data

        //- Default matching tolerance
        static const scalar defaultMatchTol_;

        //- Local matching tolerance
        const scalar matchTolerance_;

        //- Type of transformation
        transformType transform_;

        //- Offset (distance) vector from one side of the couple to the other
        mutable vectorField separation_;

        //- Face transformation tensor
        mutable tensorField forwardT_;

        //- Neighbour-cell transformation tensor
        mutable tensorField reverseT_;

        //- Are faces collocated. Either size 0,1 or length of patch.
        mutable boolList collocated_;


protected:

    // Protected Member Functions

        //- Calculate the transformation tensors
        //  smallDist : matching distance per face
        //  absTol    : absolute error in normal
        //  if transformType = unknown it first tries rotational, then
        //  translational transform
        void calcTransformTensors
        (
            const vectorField& Cf,
            const vectorField& Cr,
            const vectorField& nf,
            const vectorField& nr,
            const scalarField& smallDist,
            const scalar absTol,
            const transformType = UNKNOWN
        ) const;

        //- Initialise the calculation of the patch geometry
        virtual void initGeometry(PstreamBuffers&) = 0;

        //- Calculate the patch geometry
        virtual void calcGeometry(PstreamBuffers&) = 0;

        //- Initialise the patches for moving points
        virtual void initMovePoints(PstreamBuffers&, const pointField&) = 0;

        //- Correct patches after moving points
        virtual void movePoints(PstreamBuffers&, const pointField&) = 0;

        //- Initialise the update of the patch topology
        virtual void initUpdateMesh(PstreamBuffers&) = 0;

        //- Update of the patch topology
        virtual void updateMesh(PstreamBuffers&) = 0;


        //- Write point in OBJ format
        static void writeOBJ(Ostream& os, const point& pt);

        //- Write selected points in OBJ format
        static void writeOBJ(Ostream&, const pointField&, const labelList&);

        //- Write patch
        static void writeOBJ
        (
            const fileName&,
            const UList<face>&,
            const pointField&
        );

        //- Write edge in OBJ format
        static void writeOBJ
        (
            Ostream& os,
            const point& p0,
            const point& p1,
            label& vertI
        );

        //- Get a unique anchor point for all faces
        static pointField getAnchorPoints
        (
            const UList<face>&,
            const pointField&,
            const transformType
        );

        //- Get the number of vertices face f needs to be rotated such that
        //  its f[0] point aligns with given anchor (within tol).
        static label getRotation
        (
            const pointField& points,
            const face& f,
            const point& anchor,
            const scalar tol
        );


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
            const word& patchType,
            const transformType transform
        );

        //- Construct from dictionary
        coupledPolyPatch
        (
            const word& name,
            const dictionary& dict,
            const label index,
            const polyBoundaryMesh& bm,
            const word& patchType,
            const transformType defaultTransform = UNKNOWN
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

        //- Construct given the original patch and a map
        coupledPolyPatch
        (
            const coupledPolyPatch& pp,
            const polyBoundaryMesh& bm,
            const label index,
            const labelUList& mapAddressing,
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

            //- Type of transform
            virtual transformType transform() const
            {
                return transform_;
            }

            //- Type of transform
            //  This is currently only for use when collapsing generated
            //  meshes that can have zero area faces.
            virtual transformType& transform()
            {
                return transform_;
            }

            //- Transform a patch-based position from other side to this side
            virtual void transformPosition(pointField&) const = 0;

            //- Transform a patch-based position from other side to this side
            virtual void transformPosition(point&, const label facei) const = 0;

            //- Are the planes separated.
            virtual bool separated() const
            {
                return separation_.size();
            }

            //- If the planes are separated the separation vector.
            virtual const vectorField& separation() const
            {
                return separation_;
            }

            //- Are the cyclic planes parallel.
            virtual bool parallel() const
            {
                return forwardT_.empty();
            }

            //- Return face transformation tensor.
            virtual const tensorField& forwardT() const
            {
                return forwardT_;
            }

            //- Return neighbour-cell transformation tensor.
            virtual const tensorField& reverseT() const
            {
                return reverseT_;
            }

            //- Are faces collocated. Either size 0,1 or length of patch
            virtual const boolList& collocated() const
            {
                return collocated_;
            }

            scalar matchTolerance() const
            {
                return matchTolerance_;
            }


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
        ) = 0;

        //- Initialize ordering for primitivePatch. Does not
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
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
