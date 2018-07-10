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

Class
    Foam::oldCyclicPolyPatch

Description
    'old' style cyclic polyPatch with all faces in single patch. Does ordering
    but cannot be used to run. Writes 'type cyclic' so foamUpgradeCyclics
    can be run afterwards.
    Used to get cyclics from mesh converters that assume cyclics in single
    patch (e.g. fluent3DMeshToFoam)

SourceFiles
    oldCyclicPolyPatch.C

\*---------------------------------------------------------------------------*/

#ifndef oldCyclicPolyPatch_H
#define oldCyclicPolyPatch_H

#include "coupledPolyPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class oldCyclicPolyPatch Declaration
\*---------------------------------------------------------------------------*/

class oldCyclicPolyPatch
:
    public coupledPolyPatch
{
    // Private data

        //- Morph:angle between normals of neighbouring faces.
        //  Used to split cyclic into halves.
        scalar featureCos_;


        // For rotation

            //- Axis of rotation for rotational cyclics
            vector rotationAxis_;

            //- Point on axis of rotation for rotational cyclics
            point rotationCentre_;

        // For translation

            //- Translation vector
            vector separationVector_;


    // Private member functions

        //- Find amongst selected faces the one with the largest area
        static label findMaxArea(const pointField&, const faceList&);

        void calcTransforms();

        //- Calculate face centres
        static pointField calcFaceCentres
        (
            const UList<face>&,
            const pointField&
        );

        //- Get f[0] for all faces
        static pointField getAnchorPoints
        (
            const UList<face>&,
            const pointField&
        );

        // Face ordering

            //- Find the two parts of the faces of pp using feature edges.
            //  Returns true if successful.
            bool getGeometricHalves
            (
                const primitivePatch&,
                labelList&,
                labelList&
            ) const;

            //- Calculate geometric factors of the two halves.
            void getCentresAndAnchors
            (
                const primitivePatch&,
                const faceList& half0Faces,
                const faceList& half1Faces,

                pointField& ppPoints,
                pointField& half0Ctrs,
                pointField& half1Ctrs,
                pointField& anchors0,
                scalarField& tols
            ) const;

            //- Given matched faces matches the anchor point. Sets faceMap,
            //  rotation. Returns true if all matched.
            bool matchAnchors
            (
                const bool report,
                const primitivePatch&,
                const labelList&,
                const pointField&,
                const labelList&,
                const faceList&,
                const labelList&,
                const scalarField&,

                labelList& faceMap,
                labelList& rotation
            ) const;

            //- For rotational cases, try to find a unique face on each side
            //  of the cyclic.
            label getConsistentRotationFace
            (
                const pointField& faceCentres
            ) const;


protected:

    // Protected Member functions

        //- Initialise the calculation of the patch geometry
        virtual void initGeometry(PstreamBuffers&);

        //- Calculate the patch geometry
        virtual void calcGeometry(PstreamBuffers&);

        //- Initialise the patches for moving points
        virtual void initMovePoints(PstreamBuffers&, const pointField&);

        //- Correct patches after moving points
        virtual void movePoints(PstreamBuffers&, const pointField&);

        //- Initialise the update of the patch topology
        virtual void initUpdateMesh(PstreamBuffers&);

        //- Update of the patch topology
        virtual void updateMesh(PstreamBuffers&);

public:

    //- Runtime type information
    TypeName("oldCyclic");


    // Constructors

        //- Construct from components
        oldCyclicPolyPatch
        (
            const word& name,
            const label size,
            const label start,
            const label index,
            const polyBoundaryMesh& bm,
            const word& patchType,
            const transformType transform = UNKNOWN
        );

        //- Construct from dictionary
        oldCyclicPolyPatch
        (
            const word& name,
            const dictionary& dict,
            const label index,
            const polyBoundaryMesh& bm,
            const word& patchType
        );

        //- Construct as copy, resetting the boundary mesh
        oldCyclicPolyPatch(const oldCyclicPolyPatch&, const polyBoundaryMesh&);

        //- Construct given the original patch and resetting the
        //  face list and boundary mesh information
        oldCyclicPolyPatch
        (
            const oldCyclicPolyPatch& pp,
            const polyBoundaryMesh& bm,
            const label index,
            const label newSize,
            const label newStart
        );

        //- Construct and return a clone, resetting the boundary mesh
        virtual autoPtr<polyPatch> clone(const polyBoundaryMesh& bm) const
        {
            return autoPtr<polyPatch>(new oldCyclicPolyPatch(*this, bm));
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
                new oldCyclicPolyPatch(*this, bm, index, newSize, newStart)
            );
        }


    // Destructor

        virtual ~oldCyclicPolyPatch();


    // Member Functions

        // Access

            //- Does this side own the patch ?
            virtual bool owner() const
            {
                NotImplemented;
                return true;
            }

            //- Transform a patch-based position from other side to this side
            virtual void transformPosition(pointField& l) const
            {
                NotImplemented;
            }

            //- Transform a patch-based position from other side to this side
            virtual void transformPosition(point&, const label facei) const
            {
                NotImplemented;
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
        );

        //- Initialize ordering for primitivePatch. Does not
        //  refer to *this (except for name() and type() etc.)
        virtual void initOrder
        (
            PstreamBuffers&,
            const primitivePatch&
        ) const;

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
