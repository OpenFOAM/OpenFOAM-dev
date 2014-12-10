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
    Foam::wallBoundedParticle

Description
    Particle class that tracks on triangles of boundary faces. Use
    trackToEdge similar to trackToFace on particle.

SourceFiles
    wallBoundedParticle.C
    wallBoundedParticleTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef wallBoundedParticle_H
#define wallBoundedParticle_H

#include "particle.H"
#include "autoPtr.H"
#include "InfoProxy.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class wallBoundedParticle Declaration
\*---------------------------------------------------------------------------*/

class wallBoundedParticle
:
    public particle
{

public:

    //- Class used to pass tracking data to the trackToFace function
    template<class CloudType>
    class TrackingData
    :
        public particle::TrackingData<CloudType>
    {

    public:

        const PackedBoolList& isWallPatch_;

        // Constructors

            inline TrackingData
            (
                CloudType& cloud,
                const PackedBoolList& isWallPatch
            )
            :
                particle::TrackingData<CloudType>
                (
                    cloud
                ),
                isWallPatch_(isWallPatch)
            {}
    };


protected:

    // Protected data

        //- Particle is on mesh edge:
        //      const face& f = mesh.faces()[tetFace()]
        //      const edge e(f[meshEdgeStart_], f.nextLabel(meshEdgeStart_));
        //  Note that this real edge
        //  is also one of the edges of the face-triangle (from
        //  tetFace()+tetPt()).
        label meshEdgeStart_;

        //- Particle is on diagonal edge:
        //      const face& f = mesh.faces()[tetFace()]
        //      label faceBasePtI = mesh.tetBasePtIs()[faceI];
        //      label diagPtI = (faceBasePtI+diagEdge_)%f.size();
        //      const edge e(f[faceBasePtI], f[diagPtI]);
        label diagEdge_;


    // Protected Member Functions

        //- Construct current edge
        edge currentEdge() const;

        //- Check if inside current tet
        //void checkInside() const;

        //- Check if on current edge
        //void checkOnEdge() const;

        //- Check if point on triangle
        //void checkOnTriangle(const point&) const;

        //- Cross mesh edge into different face on same cell
        void crossEdgeConnectedFace(const edge& meshEdge);

        //- Cross diagonal edge into different triangle on same face,cell
        void crossDiagonalEdge();

        //- Track through single triangle
        scalar trackFaceTri(const vector& endPosition, label& minEdgeI);

        //- Is current triangle in the track direction
        bool isTriAlongTrack(const point& endPosition) const;


    // Patch interactions

        //- Do all patch interaction
        template<class TrackData>
        void patchInteraction(TrackData& td, const scalar trackFraction);

        //- Overridable function to handle the particle hitting a patch
        //  Executed before other patch-hitting functions
        template<class TrackData>
        bool hitPatch
        (
            const polyPatch&,
            TrackData& td,
            const label patchI,
            const scalar trackFraction,
            const tetIndices& tetIs
        );

        //- Overridable function to handle the particle hitting a wedge
        template<class TrackData>
        void hitWedgePatch
        (
            const wedgePolyPatch&,
            TrackData& td
        );

        //- Overridable function to handle the particle hitting a
        //  symmetry plane
        template<class TrackData>
        void hitSymmetryPlanePatch
        (
            const symmetryPlanePolyPatch&,
            TrackData& td
        );

        //- Overridable function to handle the particle hitting a
        //  symmetry patch
        template<class TrackData>
        void hitSymmetryPatch
        (
            const symmetryPolyPatch&,
            TrackData& td
        );

        //- Overridable function to handle the particle hitting a cyclic
        template<class TrackData>
        void hitCyclicPatch
        (
            const cyclicPolyPatch&,
            TrackData& td
        );

        //- Overridable function to handle the particle hitting a
        //- processorPatch
        template<class TrackData>
        void hitProcessorPatch
        (
            const processorPolyPatch&,
            TrackData& td
        );

        //- Overridable function to handle the particle hitting a wallPatch
        template<class TrackData>
        void hitWallPatch
        (
            const wallPolyPatch&,
            TrackData& td,
            const tetIndices&
        );

        //- Overridable function to handle the particle hitting a polyPatch
        template<class TrackData>
        void hitPatch
        (
            const polyPatch&,
            TrackData& td
        );

public:

    // Constructors

        //- Construct from components
        wallBoundedParticle
        (
            const polyMesh& c,
            const vector& position,
            const label cellI,
            const label tetFaceI,
            const label tetPtI,
            const label meshEdgeStart,
            const label diagEdge
        );

        //- Construct from Istream
        wallBoundedParticle
        (
            const polyMesh& c,
            Istream& is,
            bool readFields = true
        );

        //- Construct copy
        wallBoundedParticle(const wallBoundedParticle& p);

        //- Construct and return a clone
        autoPtr<particle> clone() const
        {
            return autoPtr<particle>(new wallBoundedParticle(*this));
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

            autoPtr<wallBoundedParticle> operator()
            (
                Istream& is
            ) const
            {
                return autoPtr<wallBoundedParticle>
                (
                    new wallBoundedParticle(mesh_, is, true)
                );
            }
        };


    // Member Functions

        // Access

            //- -1 or label of mesh edge
            inline label meshEdgeStart() const
            {
                return meshEdgeStart_;
            }

            //- -1 or diagonal edge
            inline label diagEdge() const
            {
                return diagEdge_;
            }


        // Track

            //- Equivalent of trackToFace
            template<class TrackData>
            scalar trackToEdge
            (
                TrackData& td,
                const vector& endPosition
            );


        // Info

            //- Return info proxy.
            //  Used to print particle information to a stream
            inline InfoProxy<wallBoundedParticle> info() const
            {
                return *this;
            }


        // I-O

            //- Read
            template<class CloudType>
            static void readFields(CloudType&);

            //- Write
            template<class CloudType>
            static void writeFields(const CloudType&);

            //- Write the particle data
            void write(Ostream& os, bool writeFields) const;


    // Ostream Operator

        friend Ostream& operator<<
        (
            Ostream&,
            const wallBoundedParticle&
        );

        friend Ostream& operator<<
        (
            Ostream&,
            const InfoProxy<wallBoundedParticle>&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "wallBoundedParticleTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
