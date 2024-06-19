/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2024 OpenFOAM Foundation
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
#include "barycentric.H"
#include "barycentricTensor.H"
#include "IDLList.H"
#include "pointField.H"
#include "faceList.H"
#include "OFstream.H"
#include "tetPointRef.H"
#include "FixedList.H"
#include "polyMeshTetDecomposition.H"
#include "particleMacros.H"
#include "transformer.H"

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
    // Private member data

        //- Size in bytes of the position data
        static const std::size_t sizeofPosition_;

        //- Size in bytes of the fields
        static const std::size_t sizeofFields_;


public:

    class trackingData
    {
    public:

        // Public data

            //- Reference to the mesh
            const polyMesh& mesh;

            //- Flag to indicate whether to keep particle (false = delete)
            bool keepParticle;

            //- Processor to send the particle to. -1 indicates that this
            //  particle is not to be transferred.
            label sendToProc;

            //- Patch from which to send the particle
            label sendFromPatch;

            //- Patch to which to send the particle
            label sendToPatch;

            //- Patch face to which to send the particle
            label sendToPatchFace;

            //- Position to which to send
            vector sendToPosition;

            //- Number of boundary hits that occurred during locate executions
            //  following (non-conformal) patch transfers. For reporting.
            labelList patchNLocateBoundaryHits;


        // Constructor
        template <class TrackCloudType>
        trackingData(const TrackCloudType& cloud)
        :
            mesh(cloud.pMesh()),
            keepParticle(false),
            sendToProc(-1),
            sendFromPatch(-1),
            sendToPatch(-1),
            sendToPatchFace(-1),
            sendToPosition(vector::uniform(NaN)),
            patchNLocateBoundaryHits
            (
                mesh.boundaryMesh().size()
              - mesh.globalData().processorPatches().size(),
                0
            )
        {}
    };


private:

    // Private Data

        //- Coordinates of particle
        barycentric coordinates_;

        //- Index of the cell it is in
        label celli_;

        //- Index of the face that owns the decomposed tet that the
        //  particle is in
        label tetFacei_;

        //- Index of the point on the face that defines the decomposed
        //  tet that the particle is in.  Relative to the face base
        //  point.
        label tetPti_;

        //- Face index if the particle is on a face otherwise -1
        label facei_;

        //- Fraction of time-step completed
        scalar stepFraction_;

        //- The step fraction less than the maximum reached so far. See
        //  tracking.H for details.
        scalar stepFractionBehind_;

        //- The number of tracks carried out that ended in a step fraction less
        //  than the maximum reached so far. See tracking.H for details.
        //  maxNTracksBehind_.
        label nTracksBehind_;

        //- Originating processor id
        label origProc_;

        //- Local particle id on originating processor
        label origId_;


public:

    // Static Data Members

        //- Runtime type information
        TypeName("particle");

        //- String representation of properties
        DefinePropertyList
        (
            "(coordinatesa coordinatesb coordinatesc coordinatesd) "
            "celli tetFacei tetPti facei stepFraction "
            "behind nBehind origProc origId"
        );

        //- Cumulative particle counter - used to provide unique ID
        static label particleCount_;


    // Constructors

        //- Construct from components
        particle
        (
            const polyMesh& mesh,
            const barycentric& coordinates,
            const label celli,
            const label tetFacei,
            const label tetPti,
            const label facei
        );

        //- Construct from a position and a cell, searching for the rest of the
        //  required topology
        particle
        (
            const polyMesh& mesh,
            const vector& position,
            const label celli,
            label& nLocateBoundaryHits
        );

        //- Construct from Istream
        particle(Istream&, bool readFields = true);

        //- Construct as a copy
        particle(const particle& p);

        //- Construct and return a clone
        virtual autoPtr<particle> clone() const
        {
            return autoPtr<particle>(new particle(*this));
        }

        //- Construct from Istream and return
        static autoPtr<particle> New(Istream& is)
        {
            return autoPtr<particle>(new particle(is));
        }


    //- Destructor
    virtual ~particle()
    {}


    // Member Functions

        // Access

            //- Get unique particle creation id
            inline label getNewParticleIndex() const;

            //- Return current particle coordinates
            inline const barycentric& coordinates() const;

            //- Return current cell particle is in
            inline label cell() const;

            //- Return current tet face particle is in
            inline label tetFace() const;

            //- Return current tet face particle is in
            inline label tetPt() const;

            //- Return current face particle is on otherwise -1
            inline label face() const;

            //- Return current face particle is on otherwise -1
            inline label& face();

            //- Return the fraction of time-step completed
            inline scalar stepFraction() const;

            //- Return the fraction of time-step completed
            inline scalar& stepFraction();

            //- Return the originating processor ID
            inline label origProc() const;

            //- Return the originating processor ID
            inline label& origProc();

            //- Return the particle ID on the originating processor
            inline label origId() const;

            //- Return the particle ID on the originating processor
            inline label& origId();


        // Check

            //- Return the indices of the current tet that the
            //  particle occupies.
            inline tetIndices currentTetIndices(const polyMesh& mesh) const;

            //- Return the normal of the tri on tetFacei_ for the
            //  current tet.
            inline vector normal(const polyMesh& mesh) const;

            //- Is the particle on a face?
            inline bool onFace() const;

            //- Is the particle on an internal face?
            inline bool onInternalFace(const polyMesh& mesh) const;

            //- Is the particle on a boundary face?
            inline bool onBoundaryFace(const polyMesh& mesh) const;

            //- Return the index of patch that the particle is on
            inline label patch(const polyMesh& mesh) const;

            //- Return current particle position
            inline vector position(const polyMesh& mesh) const;


        // Track

            //- Set the step fraction and clear the behind data in preparation
            //  for a new track
            inline void reset(const scalar stepFraction);

            //- Locate the particle at the given position
            bool locate
            (
                const polyMesh& mesh,
                const vector& position,
                label celli
            );

            //- Track along the displacement for a given fraction of the overall
            //  step. End when the track is complete, or when a boundary is hit.
            //  On exit, stepFraction_ will have been incremented to the current
            //  position, and facei_ will be set to the index of the boundary
            //  face that was hit, or -1 if the track completed within a cell.
            //  The proportion of the displacement still to be completed is
            //  returned.
            scalar track
            (
                const polyMesh& mesh,
                const vector& displacement,
                const scalar fraction
            );

            //- As particle::track, but stops when a new cell is reached.
            scalar trackToCell
            (
                const polyMesh& mesh,
                const vector& displacement,
                const scalar fraction
            );

            //- As particle::track, but stops when a face is hit.
            scalar trackToFace
            (
                const polyMesh& mesh,
                const vector& displacement,
                const scalar fraction
            );

            //- Hit the current face. If the current face is internal than this
            //  crosses into the next cell. If it is a boundary face then this
            //  will interact the particle with the relevant patch.
            template<class TrackCloudType>
            void hitFace
            (
                const vector& displacement,
                const scalar fraction,
                TrackCloudType& cloud,
                trackingData& td
            );

            //- Convenience function. Combines trackToFace and hitFace
            template<class TrackCloudType>
            scalar trackToAndHitFace
            (
                const vector& displacement,
                const scalar fraction,
                TrackCloudType& cloud,
                trackingData& td
            );

            //- Get the displacement from the mesh centre. Used to correct the
            //  particle position in cases with reduced dimensionality. Returns
            //  a zero vector for three-dimensional cases.
            vector deviationFromMeshCentre(const polyMesh& mesh) const;


        // Patch data

            //- Get the normal and displacement of the current patch location
            inline void patchData
            (
                const polyMesh& mesh,
                vector& normal,
                vector& displacement
            ) const;


        // Transformations

            //- Transform the physical properties of the particle
            //  according to the given transformation tensor
            virtual void transformProperties(const transformer&);


        // Transfers

            //- Make changes prior to a parallel transfer. Runs either
            //  processor or nonConformalCyclic variant below.
            template<class TrackCloudType>
            void prepareForParallelTransfer(TrackCloudType&, trackingData&);

            //- Make changes following a parallel transfer. Runs either
            //  processor or nonConformalCyclic variant below.
            template<class TrackCloudType>
            void correctAfterParallelTransfer(TrackCloudType&, trackingData&);

            //- Make changes prior to a transfer across a processor boundary.
            //  Stores the local patch face index (in facei_) so that the mesh
            //  face index can be determined on the other side.
            void prepareForProcessorTransfer(trackingData& td);

            //- Make changes following a transfer across a processor boundary.
            //  Converts the stored patch index to a mesh index. Accounts for
            //  the receiving face being reversed relative to the sending face.
            void correctAfterProcessorTransfer(trackingData& td);

            //- Make changes prior to a transfer across a non conformal cyclic
            //  boundary. Stores the receiving patch face (in facei_). Breaks
            //  the topology and stores the cartesian position.
            void prepareForNonConformalCyclicTransfer
            (
                const polyMesh& mesh,
                const label sendToPatch,
                const label sendToPatchFace,
                const vector& sendToPosition
            );

            //- Make changes following a transfer across a non conformal cyclic
            //  boundary. Locates the particle using the stored face index and
            //  cartesian position.
            void correctAfterNonConformalCyclicTransfer
            (
                const polyMesh& mesh,
                const label sendToPatch,
                labelList& patchNLocateBoundaryHits
            );


        // Patch interactions

            //- Overridable function to handle the particle hitting a patch.
            //  Executed before other patch-hitting functions.
            template<class TrackCloudType>
            bool hitPatch(TrackCloudType&, trackingData&);

            //- Overridable function to handle the particle hitting a wedgePatch
            template<class TrackCloudType>
            void hitWedgePatch(TrackCloudType&, trackingData&);

            //- Overridable function to handle the particle hitting a
            //  symmetryPlanePatch
            template<class TrackCloudType>
            void hitSymmetryPlanePatch(TrackCloudType&, trackingData&);

            //- Overridable function to handle the particle hitting a
            //  symmetryPatch
            template<class TrackCloudType>
            void hitSymmetryPatch(TrackCloudType&, trackingData&);

            //- Overridable function to handle the particle hitting a
            //  cyclicPatch
            template<class TrackCloudType>
            void hitCyclicPatch(TrackCloudType&, trackingData&);

            //- Overridable function to handle the particle hitting an
            //  nonConformalCyclicPolyPatch
            template<class TrackCloudType>
            bool hitNonConformalCyclicPatch
            (
                const vector& displacement,
                const scalar fraction,
                const label patchi,
                TrackCloudType& cloud,
                trackingData& td
            );

            //- Overridable function to handle the particle hitting a
            //  processorPatch
            template<class TrackCloudType>
            void hitProcessorPatch(TrackCloudType&, trackingData&);

            //- Overridable function to handle the particle hitting a wallPatch
            template<class TrackCloudType>
            void hitWallPatch(TrackCloudType&, trackingData&);

            //- Overridable function to handle the particle hitting a basic
            //  patch. Fall-through for the above.
            template<class TrackCloudType>
            void hitBasicPatch(TrackCloudType&, trackingData&);


        // Interaction list referral

            //- Break the topology and store the cartesian position so that the
            //  particle can be referred.
            void prepareForInteractionListReferral
            (
                const polyMesh& mesh,
                const transformer& transform
            );

            //- Correct the topology after referral. Locates the particle
            //  relative to a nearby cell/tet. The particle may end up outside
            //  this cell/tet and cannot therefore be tracked.
            void correctAfterInteractionListReferral
            (
                const polyMesh& mesh,
                const label celli
            );


        // Decompose and reconstruct

            //- Return the tet point appropriate for decomposition or
            //  reconstruction to or from the given mesh.
            label procTetPt
            (
                const polyMesh& mesh,
                const polyMesh& procMesh,
                const label procCell,
                const label procTetFace
            ) const;


        // I-O

            //- Read the fields associated with the owner cloud
            template<class TrackCloudType>
            static void readFields(TrackCloudType& c);

            //- Write the fields associated with the owner cloud
            template<class TrackCloudType>
            static void writeFields(const TrackCloudType& c);

            //- Write the particle position and cell
            void writePosition(Ostream&) const;


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
    #include "particleTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
