/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2021 OpenFOAM Foundation
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

            const PtrList<interpolation<scalar>>& vsInterp_;

            const PtrList<interpolation<vector>>& vvInterp_;

            const label UIndex_;

            bool trackForward_;

            const label nSubCycle_;

            const scalar trackLength_;

            DynamicList<vectorList>& allPositions_;

            DynamicList<scalarList>& allTimes_;

            List<DynamicList<scalarList>>& allScalars_;

            List<DynamicList<vectorList>>& allVectors_;


        // Constructors

            //- Construct from components
            trackingData
            (
                streamlinesCloud& cloud,
                const PtrList<interpolation<scalar>>& vsInterp,
                const PtrList<interpolation<vector>>& vvInterp,
                const label UIndex,
                const bool trackForward,
                const label nSubCycle,
                const scalar trackLength,
                DynamicList<List<point>>& allPositions,
                DynamicList<List<scalar>>& allTimes,
                List<DynamicList<scalarList>>& allScalars,
                List<DynamicList<vectorList>>& allVectors
            )
            :
                particle::trackingData(cloud),
                vsInterp_(vsInterp),
                vvInterp_(vvInterp),
                UIndex_(UIndex),
                trackForward_(trackForward),
                nSubCycle_(nSubCycle),
                trackLength_(trackLength),
                allPositions_(allPositions),
                allTimes_(allTimes),
                allScalars_(allScalars),
                allVectors_(allVectors)
            {}
    };


private:

    // Private Data

        //- Lifetime of particle. Particle dies when reaches 0.
        label lifeTime_;

        //- Age of the particle
        scalar age_;

        //- Sampled positions
        DynamicList<point> sampledPositions_;

        //- Sampled times
        DynamicList<scalar> sampledTimes_;

        //- Sampled scalars
        List<DynamicList<scalar>> sampledScalars_;

        //- Sampled vectors
        List<DynamicList<vector>> sampledVectors_;


    // Private Member Functions

        //- Interpolate all quantities; return interpolated velocity.
        vector interpolateFields
        (
            const trackingData&,
            const point&,
            const label celli,
            const label facei
        );


public:

    // Constructors

        //- Construct from components
        streamlinesParticle
        (
            const polyMesh& c,
            const vector& position,
            const label celli,
            const label lifeTime
        );

        //- Construct from Istream
        streamlinesParticle
        (
            const polyMesh& c,
            Istream& is,
            bool readFields = true
        );

        //- Construct copy
        streamlinesParticle(const streamlinesParticle& p);

        //- Construct and return a clone
        autoPtr<particle> clone() const
        {
            return autoPtr<particle>(new streamlinesParticle(*this));
        }

        //- Factory class to read-construct particles used for parallel transfer
        class iNew
        {
            const polyMesh& mesh_;

        public:

            iNew(const polyMesh& mesh)
            :
                mesh_(mesh)
            {}

            autoPtr<streamlinesParticle> operator()(Istream& is) const
            {
                return autoPtr<streamlinesParticle>
                (
                    new streamlinesParticle(mesh_, is, true)
                );
            }
        };


    // Member Functions

        // Tracking

            //- Track all particles to their end point
            bool move(streamlinesCloud&, trackingData&, const scalar);

            //- Overridable function to handle the particle hitting a patch
            //  Executed before other patch-hitting functions
            bool hitPatch(streamlinesCloud&, trackingData&);

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

            //- Overridable function to handle the particle hitting a
            //  cyclicAMIPatch
            void hitCyclicAMIPatch
            (
                const vector&,
                const scalar,
                streamlinesCloud&,
                trackingData&
            );

            //- Overridable function to handle the particle hitting a
            //  cyclicACMIPatch
            void hitCyclicACMIPatch
            (
                const vector&,
                const scalar,
                streamlinesCloud&,
                trackingData&
            );

            //- Overridable function to handle the particle hitting a
            //  cyclicRepeatAMIPatch
            void hitCyclicRepeatAMIPatch
            (
                const vector&,
                const scalar,
                streamlinesCloud&,
                trackingData&
            );

            //- Overridable function to handle the particle hitting a
            //- processorPatch
            void hitProcessorPatch(streamlinesCloud&, trackingData&);

            //- Overridable function to handle the particle hitting a wallPatch
            void hitWallPatch(streamlinesCloud&, trackingData&);


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
