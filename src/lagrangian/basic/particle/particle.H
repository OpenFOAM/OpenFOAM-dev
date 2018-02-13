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
#include "Cloud.H"
#include "IDLList.H"
#include "pointField.H"
#include "faceList.H"
#include "OFstream.H"
#include "tetPointRef.H"
#include "FixedList.H"
#include "polyMeshTetDecomposition.H"
#include "particleMacros.H"
#include "vectorTensorTransform.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class particle;

class polyPatch;

class cyclicPolyPatch;
class cyclicAMIPolyPatch;
class cyclicACMIPolyPatch;
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

        //- The factor by which the displacement is increased when passing
        //  through negative space. This should be slightly bigger than one.
        //  This is needed as a straight trajectory can form a closed loop
        //  through regions of overlapping positive and negative space, leading
        //  to a track which never ends. By increasing the rate of displacement
        //  through negative regions, the change in track fraction over this
        //  part of the loop no longer exactly cancels the change over the
        //  positive part, and the track therefore terminates.
        static const scalar negativeSpaceDisplacementFactor;


public:

    class trackingData
    {
    public:

        // Public data

            //- Flag to switch processor
            bool switchProcessor;

            //- Flag to indicate whether to keep particle (false = delete)
            bool keepParticle;


        // Constructor
        template <class TrackCloudType>
        trackingData(const TrackCloudType& cloud)
        {}
    };


private:

    // Private data

        //- Reference to the polyMesh database
        const polyMesh& mesh_;

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

        //- Originating processor id
        label origProc_;

        //- Local particle id on originating processor
        label origId_;


private:

    // Private Member Functions

        // Tetrahedra functions

            //- Get the vertices of the current tet
            void stationaryTetGeometry
            (
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
            barycentricTensor stationaryTetTransform() const;

            //- Get the reverse transform associated with the current tet. The
            //  conversion is detA*y = (x - centre) & T. The variables x, y and
            //  centre have the same meaning as for the forward transform. T is
            //  the transposed inverse of the forward transform tensor, A,
            //  multiplied by its determinant, detA. This separation allows
            //  the barycentric tracking algorithm to function on inverted or
            //  degenerate tetrahedra.
            void stationaryTetReverseTransform
            (
                vector& centre,
                scalar& detA,
                barycentricTensor& T
            ) const;

            //- Get the vertices of the current moving tet. Two values are
            //  returned for each vertex. The first is a constant, and the
            //  second is a linear coefficient of the track fraction.
            void movingTetGeometry
            (
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
            Pair<barycentricTensor> movingTetTransform
            (
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
            void changeTet(const label tetTriI);

            //- Change tet face within a cell. Called by changeTet.
            void changeFace(const label tetTriI);

            //- Change cell. Called when the particle hits an internal face.
            void changeCell();

            //- Put the particle on the lowest indexed patch for the current set
            //  of coincident faces. In the case of an ACMI-wall pair, this will
            //  move the particle from the wall face to the ACMI face, because
            //  ACMI patches are always listed before their associated non-
            //  overlapping patch.
            void changeToMasterPatch();


        // Geometry changes

            //- Locate the particle at the given position
            void locate
            (
                const vector& position,
                const vector* direction,
                const label celli,
                const bool boundaryFail,
                const string boundaryMsg
            );


protected:

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

        //- Overridable function to handle the particle hitting a symmetryPatch
        template<class TrackCloudType>
        void hitSymmetryPatch(TrackCloudType&, trackingData&);

        //- Overridable function to handle the particle hitting a cyclicPatch
        template<class TrackCloudType>
        void hitCyclicPatch(TrackCloudType&, trackingData&);

        //- Overridable function to handle the particle hitting a cyclicAMIPatch
        template<class TrackCloudType>
        void hitCyclicAMIPatch(TrackCloudType&, trackingData&, const vector&);

        //- Overridable function to handle the particle hitting a
        //  cyclicACMIPatch
        template<class TrackCloudType>
        void hitCyclicACMIPatch(TrackCloudType&, trackingData&, const vector&);

        //- Overridable function to handle the particle hitting a processorPatch
        template<class TrackCloudType>
        void hitProcessorPatch(TrackCloudType&, trackingData&);

        //- Overridable function to handle the particle hitting a wallPatch
        template<class TrackCloudType>
        void hitWallPatch(TrackCloudType&, trackingData&);


public:

    // Static data members

        //- Runtime type information
        TypeName("particle");

        //- String representation of properties
        DefinePropertyList
        (
            "(coordinatesa coordinatesb coordinatesc coordinatesd) "
            "celli tetFacei tetPti facei stepFraction origProc origId"
        );

        //- Cumulative particle counter - used to provode unique ID
        static label particleCount_;


    // Constructors

        //- Construct from components
        particle
        (
            const polyMesh& mesh,
            const barycentric& coordinates,
            const label celli,
            const label tetFacei,
            const label tetPti
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
        particle(const polyMesh& mesh, Istream&, bool readFields = true);

        //- Construct as a copy
        particle(const particle& p);

        //- Construct as a copy with references to a new mesh
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

            //- Return the step fraction change within the overall time-step.
            //  Returns the start value and the change as a scalar pair. Always
            //  return Pair<scalar>(0, 1), unless sub-cycling is in effect, in
            //  which case the values will reflect the span of the sub-cycle
            //  within the time-step.
            inline Pair<scalar> stepFractionSpan() const;

            //- Return the current fraction within the timestep. This differs
            //  from the stored step fraction due to sub-cycling.
            inline scalar currentTimeFraction() const;

            //- Return the indices of the current tet that the
            //  particle occupies.
            inline tetIndices currentTetIndices() const;

            //- Return the current tet transformation tensor
            inline barycentricTensor currentTetTransform() const;

            //- Return the normal of the tri on tetFacei_ for the
            //  current tet.
            inline vector normal() const;

            //- Return the normal of the tri on tetFacei_ for the
            //  current tet at the start of the timestep, i.e. based
            //  on oldPoints
            inline vector oldNormal() const;

            //- Is the particle on a face?
            inline bool onFace() const;

            //- Is the particle on an internal face?
            inline bool onInternalFace() const;

            //- Is the particle on a boundary face?
            inline bool onBoundaryFace() const;

            //- Return the index of patch that the particle is on
            inline label patch() const;

            //- Return current particle position
            inline vector position() const;


    // Track

        //- Track along the displacement for a given fraction of the overall
        //  step. End when the track is complete, or when a boundary is hit.
        //  On exit, stepFraction_ will have been incremented to the current
        //  position, and facei_ will be set to the index of the boundary
        //  face that was hit, or -1 if the track completed within a cell.
        //  The proportion of the displacement still to be completed is
        //  returned.
        scalar track
        (
            const vector& displacement,
            const scalar fraction
        );

        //- As particle::track, but also stops on internal faces.
        scalar trackToFace
        (
            const vector& displacement,
            const scalar fraction
        );

        //- As particle::trackToFace, but also stops on tet triangles. On
        //  exit, tetTriI is set to the index of the tet triangle that was
        //  hit, or -1 if the end position was reached.
        scalar trackToTri
        (
            const vector& displacement,
            const scalar fraction,
            label& tetTriI
        );

        //- As particle::trackToTri, but for stationary meshes
        scalar trackToStationaryTri
        (
            const vector& displacement,
            const scalar fraction,
            label& tetTriI
        );

        //- As particle::trackToTri, but for moving meshes
        scalar trackToMovingTri
        (
            const vector& displacement,
            const scalar fraction,
            label& tetTriI
        );

        //- Hit the current face. If the current face is internal than this
        //  crosses into the next cell. If it is a boundary face then this will
        //  interact the particle with the relevant patch.
        template<class TrackCloudType>
        void hitFace
        (
            const vector& direction,
            TrackCloudType& cloud,
            trackingData& td
        );

        //- Convenience function. Cobines trackToFace and hitFace
        template<class TrackCloudType>
        void trackToAndHitFace
        (
            const vector& direction,
            const scalar fraction,
            TrackCloudType& cloud,
            trackingData& td
        );

        //- Get the displacement from the mesh centre. Used to correct the
        //  particle position in cases with reduced dimensionality. Returns a
        //  zero vector for three-dimensional cases.
        vector deviationFromMeshCentre() const;


    // Patch data

        //- Get the normal and velocity of the current patch location
        void patchData(vector& n, vector& U) const;


    // Transformations

        //- Transform the physical properties of the particle
        //  according to the given transformation tensor
        virtual void transformProperties(const tensor& T);

        //- Transform the physical properties of the particle
        //  according to the given separation vector
        virtual void transformProperties(const vector& separation);


    // Parallel transfer

        //- Convert global addressing to the processor patch local equivalents
        void prepareForParallelTransfer();

        //- Convert processor patch addressing to the global equivalents
        //  and set the celli to the face-neighbour
        void correctAfterParallelTransfer(const label patchi, trackingData& td);


    // Interaction list referral

        //- Break the topology and store the particle position so that the
        //  particle can be referred.
        void prepareForInteractionListReferral
        (
            const vectorTensorTransform& transform
        );

        //- Correct the topology after referral. The particle may still be
        //  outside the stored tet and therefore not track-able.
        void correctAfterInteractionListReferral(const label celli);


    // Decompose and reconstruct

        //- Return the tet point approproate for decomposition or reconstruction
        //  to or from the given mesh.
        label procTetPt
        (
            const polyMesh& procMesh,
            const label procCell,
            const label procTetFace
        ) const;


    // Mapping

        //- Map after a topology change
        void autoMap(const vector& position, const mapPolyMesh& mapper);


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
