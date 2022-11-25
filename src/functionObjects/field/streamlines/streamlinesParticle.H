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
    Foam::streamlinesParticle

Description
    Particle class that samples fields as it passes through. Used in streamlines
    calculation.

SourceFiles
    streamlinesParticle.C

\*---------------------------------------------------------------------------*/

#ifndef streamlinesParticle_H
#define streamlinesParticle_H

#include "particle.H"
#include "Cloud.H"
#include "autoPtr.H"
#include "interpolation.H"
#include "vectorList.H"
#include "DynamicField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

class streamlinesParticle;
class streamlinesCloud;

Ostream& operator<<(Ostream&, const streamlinesParticle&);

/*---------------------------------------------------------------------------*\
                     Class streamlinesParticle Declaration
\*---------------------------------------------------------------------------*/

class streamlinesParticle
:
    public particle
{
public:

    class trackingData
    :
        public particle::trackingData
    {
    public:

        // Public data

            #define DeclareTypeInterpolator(Type, nullArg) \
                const PtrList<interpolation<Type>>& Type##Interp_;
            FOR_ALL_FIELD_TYPES(DeclareTypeInterpolator);
            #undef DeclareTypeInterpolator

            const interpolation<vector>& UInterp_;

            bool trackForward_;

            bool trackOutside_;

            const label nSubCycle_;

            const scalar trackLength_;

            label tracki;

            DynamicField<point>& allPositions_;

            DynamicField<label>& allTracks_;

            DynamicField<label>& allTrackParts_;

            DynamicField<scalar>& allAges_;

            #define DeclareAllTypes(Type, nullArg) \
                List<DynamicField<Type>>& all##Type##s_;
            FOR_ALL_FIELD_TYPES(DeclareAllTypes);
            #undef DeclareAllTypes


        // Constructors

            //- Construct from components
            trackingData
            (
                streamlinesCloud& cloud,
                #define TypeInterpolatorArg(Type, nullArg) \
                    const PtrList<interpolation<Type>>& Type##Interp,
                FOR_ALL_FIELD_TYPES(TypeInterpolatorArg)
                #undef TypeInterpolatorArg
                const interpolation<vector>& UInterp,
                const bool trackForward,
                const bool trackOutside,
                const label nSubCycle,
                const scalar trackLength,
                DynamicField<point>& allPositions,
                DynamicField<label>& allTracks,
                DynamicField<label>& allTrackParts,
                DynamicField<scalar>& allAges
                #define AllTypesArg(Type, nullArg) \
                    , List<DynamicField<Type>>& all##Type##s
                FOR_ALL_FIELD_TYPES(AllTypesArg)
                #undef AllTypesArg
            )
            :
                particle::trackingData(cloud),
                #define TypeInterpolatorInit(Type, nullArg) \
                    Type##Interp_(Type##Interp),
                FOR_ALL_FIELD_TYPES(TypeInterpolatorInit)
                #undef TypeInterpolatorInit
                UInterp_(UInterp),
                trackForward_(trackForward),
                trackOutside_(trackOutside),
                nSubCycle_(nSubCycle),
                trackLength_(trackLength),
                allPositions_(allPositions),
                allTracks_(allTracks),
                allTrackParts_(allTrackParts),
                allAges_(allAges)
                #define AllTypesInit(Type, nullArg) \
                    , all##Type##s_(all##Type##s)
                FOR_ALL_FIELD_TYPES(AllTypesInit)
                #undef AllTypesInit
            {}
    };


private:

    // Private Data

        //- Lifetime of particle. Particle dies when reaches 0.
        label lifeTime_;

        //- Index of the track
        label trackIndex_;

        //- Index of the part of the track
        label trackPartIndex_;

        //- Age of the particle
        scalar age_;

        //- Current compound transform
        transformer transform_;

        //- Sampled positions
        DynamicField<point> sampledPositions_;

        //- Sampled ages
        DynamicField<scalar> sampledAges_;

        //- Sampled types
        #define DeclareSampledTypes(Type, nullArg) \
            List<DynamicField<Type>> sampled##Type##s_;
        FOR_ALL_FIELD_TYPES(DeclareSampledTypes);
        #undef DeclareSampledTypes


    // Private Member Functions

        //- Interpolate all quantities; return interpolated velocity.
        vector interpolateFields
        (
            const trackingData&,
            const point&,
            const label celli,
            const label facei
        );

        //- End the current track
        void endTrack(trackingData&);


public:

    // Constructors

        //- Construct from components
        streamlinesParticle
        (
            const polyMesh& mesh,
            const vector& position,
            const label celli,
            const label lifeTime,
            const label trackIndex
        );

        //- Construct from Istream
        streamlinesParticle(Istream& is, bool readFields = true);

        //- Construct copy
        streamlinesParticle(const streamlinesParticle& p);

        //- Construct and return a clone
        autoPtr<particle> clone() const
        {
            return autoPtr<particle>(new streamlinesParticle(*this));
        }

        //- Construct from Istream and return
        static autoPtr<streamlinesParticle> New(Istream& is)
        {
            return autoPtr<streamlinesParticle>(new streamlinesParticle(is));
        }


    // Member Functions

        // Tracking

            //- Track all particles to their end point
            bool move(streamlinesCloud&, trackingData&);

            //- Overridable function to handle the particle hitting a wedge
            void hitWedgePatch(streamlinesCloud&, trackingData&);

            //- Overridable function to handle the particle hitting a
            //  symmetry plane
            void hitSymmetryPlanePatch(streamlinesCloud&, trackingData&);

            //- Overridable function to handle the particle hitting a
            //  symmetry patch
            void hitSymmetryPatch(streamlinesCloud&, trackingData&);

            //- Overridable function to handle the particle hitting a cyclic
            void hitCyclicPatch(streamlinesCloud&, trackingData&);

            //- Overridable function to handle the particle hitting an
            //  nonConformalCyclicPolyPatch
            bool hitNonConformalCyclicPatch
            (
                const vector& displacement,
                const scalar fraction,
                const label patchi,
                streamlinesCloud& cloud,
                trackingData& td
            );

            //- Overridable function to handle the particle hitting a
            //  processorPatch
            void hitProcessorPatch(streamlinesCloud&, trackingData&);

            //- Overridable function to handle the particle hitting a wallPatch
            void hitWallPatch(streamlinesCloud&, trackingData&);


        // Transformations

            //- Transform the physical properties of the particle
            //  according to the given transformation tensor
            virtual void transformProperties(const transformer&);


        // I-O

            //- Read
            static void readFields(Cloud<streamlinesParticle>&);

            //- Write
            static void writeFields(const Cloud<streamlinesParticle>&);


    // Ostream Operator

        friend Ostream& operator<<(Ostream&, const streamlinesParticle&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
