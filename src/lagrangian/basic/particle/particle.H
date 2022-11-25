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

        //- The value of nTracksBehind_ at which tracking is abandoned. See the
        //  description of nTracksBehind_.
        static const label maxNTracksBehind_;


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


        // Constructor
        template <class TrackCloudType>
        trackingData(const TrackCloudType& cloud)
        :
            mesh(cloud.pMesh()),
            keepParticle(false),
            sendToProc(-1),
            sendFromPatch(-1),
            sendToPatch(-1),
            sendToPatchFace(-1)
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

        //- The step fraction less than the maximum reached so far
        scalar stepFractionBehind_;

        //- The number of tracks carried out that ended in a step fraction less
        //  than the maximum reached so far. Once this reaches
        //  maxNTracksBehind_, tracking is abandoned for the current step.
        //
        //  This is needed because when tetrahedra are inverted a straight
        //  trajectory can form a closed loop through regions of overlapping
        //  positive and negative space. Without this break clause, such loops
        //  can result in a valid track which never ends.
        //
        //  Because the test is susceptible to round off error, a track of zero
        //  length will also increment the counter. As such, it is important
        //  that maxNTracksBehind_ is set large enough so that valid small
        //  tracks do not result in the track being abandoned. The largest
        //  number of valid small tracks that are likely to be performed
        //  sequentially is equal to the number of tetrahedra that can meet at
        //  a point. An estimate of this number is therefore used to set
        //  maxNTracksBehind_.
        label nTracksBehind_;

        //- Originating processor id
        label origProc_;

        //- Local particle id on originating processor
        label origId_;


private:

    // Private Member Functions

        // Tetrahedra functions

            //- Get the vertices of the current tet
            inline void stationaryTetGeometry
            (
                const polyMesh& mesh,
                vector& centre,
                vector& base,
                vector& vertex1,
                vector& vertex2
            ) const;

            //- Get the transformation associated with the current tet. This
            //  will convert a barycentric position within the tet to a
            //  cartesian position in the global coordinate system. The
            //  conversion is x = A & y, where x is the cartesian position, y is
            //  the barycentric position and A is the transformation tensor.
            inline barycentricTensor stationaryTetTransform
            (
                const polyMesh& mesh
            ) const;

            //- Get the reverse transform associated with the current tet. The
            //  conversion is detA*y = (x - centre) & T. The variables x, y and
            //  centre have the same meaning as for the forward transform. T is
            //  the transposed inverse of the forward transform tensor, A,
            //  multiplied by its determinant, detA. This separation allows
            //  the barycentric tracking algorithm to function on inverted or
            //  degenerate tetrahedra.
            void stationaryTetReverseTransform
            (
                const polyMesh& mesh,
                vector& centre,
                scalar& detA,
                barycentricTensor& T
            ) const;

            //- Get the vertices of the current moving tet. Two values are
            //  returned for each vertex. The first is a constant, and the
            //  second is a linear coefficient of the track fraction.
            inline void movingTetGeometry
            (
                const polyMesh& mesh,
                const scalar endStepFraction,
                Pair<vector>& centre,
                Pair<vector>& base,
                Pair<vector>& vertex1,
                Pair<vector>& vertex2
            ) const;

            //- Get the transformation associated with the current, moving, tet.
            //  This is of the same form as for the static case. As with the
            //  moving geometry, a linear function of the tracking fraction is
            //  returned for each component.
            inline Pair<barycentricTensor> movingTetTransform
            (
                const polyMesh& mesh,
                const scalar endStepFraction
            ) const;

            //- Get the reverse transformation associated with the current,
            //  moving, tet. This is of the same form as for the static case. As
            //  with the moving geometry, a function of the tracking fraction is
            //  returned for each component. The functions are higher order than
            //  for the forward transform; the determinant is cubic, and the
            //  tensor is quadratic.
            void movingTetReverseTransform
            (
                const polyMesh& mesh,
                const scalar endStepFraction,
                Pair<vector>& centre,
                FixedList<scalar, 4>& detA,
                FixedList<barycentricTensor, 3>& T
            ) const;


        // Transformations

            //- Reflection transform. Corrects the coordinates when the particle
            //  moves between two tets which share a base vertex, but for which
            //  the other two non cell-centre vertices are reversed. All hits
            //  which retain the same face behave this way, as do face hits.
            void reflect();

            //- Rotation transform. Corrects the coordinates when the particle
            //  moves between two tets with different base vertices, but are
            //  otherwise similarly oriented. Hits which change the face within
            //  the cell make use of both this and the reflect transform.
            void rotate(const bool direction);


        // Topology changes

            //- Change tet within a cell. Called after a triangle is hit.
            void changeTet(const polyMesh& mesh, const label tetTriI);

            //- Change tet face within a cell. Called by changeTet.
            void changeFace(const polyMesh& mesh, const label tetTriI);

            //- Change cell. Called when the particle hits an internal face.
            void changeCell(const polyMesh& mesh);


        // Geometry changes

            //- Locate the particle at the given position
            void locate
            (
                const polyMesh& mesh,
                const vector& position,
                label celli,
                const bool boundaryFail,
                const string boundaryMsg
            );


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
            const label celli
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
            inline label getNewParticleID() const;

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

            //- Return the current tet transformation tensor
            inline barycentricTensor currentTetTransform
            (
                const polyMesh& mesh
            ) const;

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

            //- As particle::trackToFace, but stops when a tet triangle is hit.
            //  On exit, tetTriI is set to the index of the tet triangle that
            //  was hit, or -1 if the end position was reached.
            scalar trackToTri
            (
                const polyMesh& mesh,
                const vector& displacement,
                const scalar fraction,
                label& tetTriI
            );

            //- As particle::trackToTri, but for stationary meshes
            scalar trackToStationaryTri
            (
                const polyMesh& mesh,
                const vector& displacement,
                const scalar fraction,
                label& tetTriI
            );

            //- As particle::trackToTri, but for moving meshes
            scalar trackToMovingTri
            (
                const polyMesh& mesh,
                const vector& displacement,
                const scalar fraction,
                label& tetTriI
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
                const label sendToPatchFace
            );

            //- Make changes following a transfer across a non conformal cyclic
            //  boundary. Locates the particle using the stored face index and
            //  cartesian position.
            void correctAfterNonConformalCyclicTransfer
            (
                const polyMesh& mesh,
                const label sendToPatch
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


        // Mapping

            //- Map after a mesh change
            void map
            (
                const polyMesh& mesh,
                const point& position,
                const label celli
            );


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
