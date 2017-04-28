/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2017 OpenFOAM Foundation
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
#include "tetrahedron.H"
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

            //- Get indices into the current face for the face-bound vertices of
            //  the current tet.
            void tetFaceIndices
            (
                label& baseI,
                label& vertex1I,
                label& vertex2I
            ) const;

            //- Get indices into the mesh points for the face-bound vertices of
            //  the current tet.
            void tetMeshIndices
            (
                label& baseI,
                label& vertex1I,
                label& vertex2I
            ) const;

            //- Get the vertices of the current tet
            void tetGeometry
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
            barycentricTensor tetTransform() const;

            //- Get the reverse transform associated with the current tet. The
            //  conversion is detA*y = (x - centre) & T. The variables x, y and
            //  centre have the same meaning as for the forward transform. T is
            //  the transposed inverse of the forward transform tensor, A,
            //  multiplied by its determinant, detA. This separation allows
            //  the barycentric tracking algorithm to function on inverted or
            //  degenerate tetrahedra.
            void tetReverseTransform
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
            const label patchi,
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
        DefinePropertyList
        (
            "(Px Py Pz) celli tetFacei tetPti "
            "facei stepFraction origProc origId"
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

            //- Return current particle coordinates
            inline const barycentric& coordinates() const;

            //- Return current cell particle is in
            inline label cell() const;

            //- Return current tet face particle is in
            inline label tetFace() const;

            //- Return current tet face particle is in
            inline label tetPt() const;

            //- Return the indices of the current tet that the
            //  particle occupies.
            inline tetIndices currentTetIndices() const;

            //- Return the geometry of the current tet that the
            //  particle occupies.
            inline tetPointRef currentTet() const;

            //- Return the normal of the tri on tetFacei_ for the
            //  current tet.
            inline vector normal() const;

            //- Return the normal of the tri on tetFacei_ for the
            //  current tet at the start of the timestep, i.e. based
            //  on oldPoints
            inline vector oldNormal() const;

            //- Return current face particle is on otherwise -1
            inline label face() const;

            //- Return the particle current time
            inline scalar currentTime() const;


        // Check

            //- Is the particle on a face?
            inline bool onFace() const;

            //- Is the particle on an internal face?
            inline bool onInternalFace() const;

            //- Is the particle on a boundary face?
            inline bool onBoundaryFace() const;

            //- Which patch is particle on // <-- !!!
            inline label patch(const label facei) const;

            //- Which face of this patch is this particle on // <-- !!!
            inline label patchFace(const label patchi, const label facei) const;

            //- Return current particle position
            inline vector position() const;

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

            //- As non-templated particle::trackToFace, but with additional
            //  boundary handling.
            template<class TrackData>
            void trackToFace
            (
                const vector& displacement,
                const scalar fraction,
                TrackData& td
            );

            //- Set the constrained components of the particle position to the
            //  mesh centre.
            void constrainToMeshCentre();


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
        void prepareForParallelTransfer(const label patchi, TrackData& td);

        //- Convert processor patch addressing to the global equivalents
        //  and set the celli to the face-neighbour
        template<class TrackData>
        void correctAfterParallelTransfer(const label patchi, TrackData& td);


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
        template<class CloudType>
        static void readFields(CloudType& c);

        //- Write the fields associated with the owner cloud
        template<class CloudType>
        static void writeFields(const CloudType& c);

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
