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
    Foam::particle

Description
    Base particle class

\*---------------------------------------------------------------------------*/

#ifndef particle_H
#define particle_H

#include "vector.H"
#include "Cloud.H"
#include "IDLList.H"
#include "pointField.H"
#include "faceList.H"
#include "OFstream.H"
#include "tetrahedron.H"
#include "FixedList.H"
#include "polyMeshTetDecomposition.H"
#include "particleMacros.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class particle;

class polyPatch;

class cyclicPolyPatch;
class processorPolyPatch;
class symmetryPlanePolyPatch;
class symmetryPolyPatch;
class wallPolyPatch;
class wedgePolyPatch;

// Forward declaration of friend functions and operators

Ostream& operator<<
(
    Ostream&,
    const particle&
);

bool operator==(const particle&, const particle&);

bool operator!=(const particle&, const particle&);

/*---------------------------------------------------------------------------*\
                          Class Particle Declaration
\*---------------------------------------------------------------------------*/

class particle
:
    public IDLList<particle>::link
{
public:

    template<class CloudType>
    class TrackingData
    {
        // Private data

            //- Reference to the cloud containing (this) particle
            CloudType& cloud_;


    public:

        // Public data

            typedef CloudType cloudType;

            //- Flag to switch processor
            bool switchProcessor;

            //- Flag to indicate whether to keep particle (false = delete)
            bool keepParticle;


        // Constructor
        TrackingData(CloudType& cloud)
        :
            cloud_(cloud)
        {}


        // Member functions

            //- Return a reference to the cloud
            CloudType& cloud()
            {
                return cloud_;
            }
    };


protected:

    // Protected data

        //- Reference to the polyMesh database
        const polyMesh& mesh_;

        //- Position of particle
        vector position_;

        //- Index of the cell it is in
        label cellI_;

        //- Face index if the particle is on a face otherwise -1
        label faceI_;

        //- Fraction of time-step completed
        scalar stepFraction_;

        //- Index of the face that owns the decomposed tet that the
        //  particle is in
        label tetFaceI_;

        //- Index of the point on the face that defines the decomposed
        //  tet that the particle is in.  Relative to the face base
        //  point.
        label tetPtI_;

        //- Originating processor id
        label origProc_;

        //- Local particle id on originating processor
        label origId_;


    // Private Member Functions

        //- Find the tet tri faces between position and tet centre
        void findTris
        (
            const vector& position,
            DynamicList<label>& faceList,
            const tetPointRef& tet,
            const FixedList<vector, 4>& tetAreas,
            const FixedList<label, 4>& tetPlaneBasePtIs,
            const scalar tol
        ) const;

        //- Find the lambda value for the line to-from across the
        //  given tri face, where p = from + lambda*(to - from)
        inline scalar tetLambda
        (
            const vector& from,
            const vector& to,
            const label triI,
            const vector& tetArea,
            const label tetPlaneBasePtI,
            const label cellI,
            const label tetFaceI,
            const label tetPtI,
            const scalar tol
        ) const;

        //- Find the lambda value for a moving tri face
        inline scalar movingTetLambda
        (
            const vector& from,
            const vector& to,
            const label triI,
            const vector& tetArea,
            const label tetPlaneBasePtI,
            const label cellI,
            const label tetFaceI,
            const label tetPtI,
            const scalar tol
        ) const;

        //- Modify the tet owner data by crossing triI
        inline void tetNeighbour(label triI);

        //- Cross the from the given face across the given edge of the
        //  given cell to find the resulting face and tetPtI
        inline void crossEdgeConnectedFace
        (
            const label& cellI,
            label& tetFaceI,
            label& tetPtI,
            const edge& e
        );

        //- Hit wall faces in the current cell if the
        //- wallImpactDistance is non-zero.  They may not be in
        //- different tets to the current.
        template<class CloudType>
        inline void hitWallFaces
        (
            const CloudType& td,
            const vector& from,
            const vector& to,
            scalar& lambdaMin,
            tetIndices& closestTetIs
        );


    // Patch interactions

        //- Overridable function to handle the particle hitting a face
        template<class TrackData>
        void hitFace(TrackData& td);

        //- Overridable function to handle the particle hitting a
        //  patch.  Executed before other patch-hitting functions.
        //  trackFraction is passed in to allow mesh motion to
        //  interpolate in time to the correct face state.
        template<class TrackData>
        bool hitPatch
        (
            const polyPatch&,
            TrackData& td,
            const label patchI,
            const scalar trackFraction,
            const tetIndices& tetIs
        );

        //- Overridable function to handle the particle hitting a wedgePatch
        template<class TrackData>
        void hitWedgePatch(const wedgePolyPatch&, TrackData& td);

        //- Overridable function to handle the particle hitting a
        //  symmetryPlanePatch
        template<class TrackData>
        void hitSymmetryPlanePatch
        (
            const symmetryPlanePolyPatch&,
            TrackData& td
        );

        //- Overridable function to handle the particle hitting a
        //  symmetryPatch
        template<class TrackData>
        void hitSymmetryPatch(const symmetryPolyPatch&, TrackData& td);

        //- Overridable function to handle the particle hitting a cyclicPatch
        template<class TrackData>
        void hitCyclicPatch(const cyclicPolyPatch&, TrackData& td);

        //- Overridable function to handle the particle hitting a cyclicAMIPatch
        template<class TrackData>
        void hitCyclicAMIPatch
        (
            const cyclicAMIPolyPatch&,
            TrackData& td,
            const vector& direction
        );

        //- Overridable function to handle the particle hitting a
        //  processorPatch
        template<class TrackData>
        void hitProcessorPatch(const processorPolyPatch&, TrackData& td);

        //- Overridable function to handle the particle hitting a wallPatch
        template<class TrackData>
        void hitWallPatch
        (
            const wallPolyPatch&,
            TrackData& td,
            const tetIndices& tetIs
        );

        //- Overridable function to handle the particle hitting a
        //  general patch
        template<class TrackData>
        void hitPatch(const polyPatch&, TrackData& td);


public:

    // Static data members

        //- Runtime type information
        TypeName("particle");

        //- String representation of properties
        DefinePropertyList("(Px Py Pz) cellI tetFaceI tetPtI origProc origId");

        //- Cumulative particle counter - used to provode unique ID
        static label particleCount_;

        //- Fraction of distance to tet centre to move a particle to
        // 'rescue' it from a tracking problem
        static const scalar trackingCorrectionTol;

        //- Fraction of the cell volume to use in determining tolerance values
        //  for the denominator and numerator of lambda
        static const scalar lambdaDistanceToleranceCoeff;

        //- Minimum stepFraction tolerance
        static const scalar minStepFractionTol;


    // Constructors

        //- Construct from components
        particle
        (
            const polyMesh& mesh,
            const vector& position,
            const label cellI,
            const label tetFaceI,
            const label tetPtI
        );

        //- Construct from components, tetFaceI_ and tetPtI_ are not
        //  supplied so they will be deduced by a search
        particle
        (
            const polyMesh& mesh,
            const vector& position,
            const label cellI,
            bool doCellFacePt = true
        );

        //- Construct from Istream
        particle(const polyMesh& mesh, Istream&, bool readFields = true);

        //- Construct as a copy
        particle(const particle& p);

        //- Construct as a copy with refernce to a new mesh
        particle(const particle& p, const polyMesh& mesh);

        //- Construct a clone
        virtual autoPtr<particle> clone() const
        {
            return autoPtr<particle>(new particle(*this));
        }


        //- Factory class to read-construct particles used for
        //  parallel transfer
        class iNew
        {
            const polyMesh& mesh_;

        public:

            iNew(const polyMesh& mesh)
            :
                mesh_(mesh)
            {}

            autoPtr<particle> operator()(Istream& is) const
            {
                return autoPtr<particle>(new particle(mesh_, is, true));
            }
        };


    //- Destructor
    virtual ~particle()
    {}


    // Member Functions

        // Access

            //- Get unique particle creation id
            inline label getNewParticleID() const;

            //- Return the mesh database
            inline const polyMesh& mesh() const;

            //- Return current particle position
            inline const vector& position() const;

            //- Return current particle position
            inline vector& position();

            //- Return current cell particle is in
            inline label& cell();

            //- Return current cell particle is in
            inline label cell() const;

            //- Return current tet face particle is in
            inline label& tetFace();

            //- Return current tet face particle is in
            inline label tetFace() const;

            //- Return current tet face particle is in
            inline label& tetPt();

            //- Return current tet face particle is in
            inline label tetPt() const;

            //- Return the indices of the current tet that the
            //  particle occupies.
            inline tetIndices currentTetIndices() const;

            //- Return the geometry of the current tet that the
            //  particle occupies.
            inline tetPointRef currentTet() const;

            //- Return the normal of the tri on tetFaceI_ for the
            //  current tet.
            inline vector normal() const;

            //- Return the normal of the tri on tetFaceI_ for the
            //  current tet at the start of the timestep, i.e. based
            //  on oldPoints
            inline vector oldNormal() const;

            //- Return current face particle is on otherwise -1
            inline label& face();

            //- Return current face particle is on otherwise -1
            inline label face() const;

            //- Return the impact model to be used, soft or hard (default).
            inline bool softImpact() const;

            //- Return the particle current time
            inline scalar currentTime() const;


        // Check

            //- Check the stored cell value (setting if necessary) and
            //  initialise the tetFace and tetPt values
            inline void initCellFacePt();

            //- Is the particle on the boundary/(or outside the domain)?
            inline bool onBoundary() const;

            //- Is this global face an internal face?
            inline bool internalFace(const label faceI) const;

            //- Is this global face a boundary face?
            inline bool boundaryFace(const label faceI) const;

            //- Which patch is particle on
            inline label patch(const label faceI) const;

            //- Which face of this patch is this particle on
            inline label patchFace
            (
                const label patchI,
                const label faceI
            ) const;

            //- Return the fraction of time-step completed
            inline scalar& stepFraction();

            //-  Return the fraction of time-step completed
            inline scalar stepFraction() const;

            //- Return const access to the originating processor id
            inline label origProc() const;

            //- Return the originating processor id for manipulation
            inline label& origProc();

            //- Return const access to  the particle id on originating processor
            inline label origId() const;

            //- Return the particle id on originating processor for manipulation
            inline label& origId();


        // Track

            //- Track particle to end of trajectory
            //  or until it hits the boundary.
            //  On entry 'stepFraction()' should be set to the fraction of the
            //  time-step at which the tracking starts and on exit it contains
            //  the fraction of the time-step completed.
            //  Returns the boundary face index if the track stops at the
            //  boundary, -1 otherwise.
            template<class TrackData>
            label track(const vector& endPosition, TrackData& td);

            //- Track particle to a given position and returns 1.0 if the
            //  trajectory is completed without hitting a face otherwise
            //  stops at the face and returns the fraction of the trajectory
            //  completed.
            //  on entry 'stepFraction()' should be set to the fraction of the
            //  time-step at which the tracking starts.
            template<class TrackData>
            scalar trackToFace(const vector& endPosition, TrackData& td);

            //- Return the index of the face to be used in the interpolation
            //  routine
            inline label faceInterpolation() const;


    // Transformations

        //- Transform the physical properties of the particle
        //  according to the given transformation tensor
        virtual void transformProperties(const tensor& T);

        //- Transform the physical properties of the particle
        //  according to the given separation vector
        virtual void transformProperties(const vector& separation);

        //- The nearest distance to a wall that
        //  the particle can be in the n direction
        virtual scalar wallImpactDistance(const vector& n) const;


    // Parallel transfer

        //- Convert global addressing to the processor patch
        //  local equivalents
        template<class TrackData>
        void prepareForParallelTransfer(const label patchI, TrackData& td);

        //- Convert processor patch addressing to the global equivalents
        //  and set the cellI to the face-neighbour
        template<class TrackData>
        void correctAfterParallelTransfer(const label patchI, TrackData& td);


    // I-O

        //- Read the fields associated with the owner cloud
        template<class CloudType>
        static void readFields(CloudType& c);

        //- Write the fields associated with the owner cloud
        template<class CloudType>
        static void writeFields(const CloudType& c);

        //- Write the particle data
        void write(Ostream& os, bool writeFields) const;


    // Friend Operators

        friend Ostream& operator<<(Ostream&, const particle&);

        friend bool operator==(const particle& pA, const particle& pB);

        friend bool operator!=(const particle& pA, const particle& pB);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "particleI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "particleTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
