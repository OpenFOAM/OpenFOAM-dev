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
    Foam::trackedParticle

Description
    Particle class that marks cells it passes through. Used to mark cells
    visited by feature edges.

SourceFiles
    trackedParticle.C

\*---------------------------------------------------------------------------*/

#ifndef trackedParticle_H
#define trackedParticle_H

#include "particle.H"
#include "autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class trackedParticleCloud;


// Forward declaration of friend functions and operators

class trackedParticle;

Ostream& operator<<(Ostream&, const trackedParticle&);


/*---------------------------------------------------------------------------*\
                     Class trackedParticle Declaration
\*---------------------------------------------------------------------------*/

class trackedParticle
:
    public particle
{
    // Private Data

        //- Start point to track from
        point start_;

        //- End point to track to
        point end_;

        //- Level of this particle
        label level_;

        //- Passive label (used to store feature edge mesh)
        label i_;

        //- Passive label (used to store feature edge point)
        label j_;

        //- Passive label (used to store feature edge label)
        label k_;


public:

    friend class Cloud<trackedParticle>;

    //- Class used to pass tracking data to the trackToFace function
    class trackingData
    :
        public particle::trackingData
    {
    public:

        labelList& maxLevel_;

        List<PackedBoolList>& featureEdgeVisited_;


        // Constructors

            trackingData
            (
                Cloud<trackedParticle>& cloud,
                labelList& maxLevel,
                List<PackedBoolList>& featureEdgeVisited
            )
            :
                particle::trackingData(cloud),
                maxLevel_(maxLevel),
                featureEdgeVisited_(featureEdgeVisited)
            {}

    };



    // Constructors

        //- Construct from components
        trackedParticle
        (
            const polyMesh& mesh,
            const barycentric& coordinates,
            const label celli,
            const label tetFacei,
            const label tetPtI,
            const point& end,
            const label level,
            const label i,
            const label j,
            const label k
        );

        //- Construct from a position and a cell, searching for the rest of the
        //  required topology
        trackedParticle
        (
            const polyMesh& mesh,
            const vector& position,
            const label celli,
            const point& end,
            const label level,
            const label i,
            const label j,
            const label k
        );

        //- Construct from Istream
        trackedParticle
        (
            const polyMesh& mesh,
            Istream& is,
            bool readFields = true
        );

        //- Construct and return a clone
        autoPtr<particle> clone() const
        {
            return autoPtr<particle>(new trackedParticle(*this));
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

            autoPtr<trackedParticle> operator()(Istream& is) const
            {
                return autoPtr<trackedParticle>
                (
                    new trackedParticle(mesh_, is, true)
                );
            }
        };


    // Member Functions

        //- Point to track from
        point& start()
        {
            return start_;
        }

        //- Point to track to
        point& end()
        {
            return end_;
        }

        //- Transported label
        label i() const
        {
            return i_;
        }

        //- Transported label
        label& i()
        {
            return i_;
        }

        //- Transported label
        label j() const
        {
            return j_;
        }

        //- Transported label
        label& j()
        {
            return j_;
        }

        //- Transported label
        label k() const
        {
            return k_;
        }

        //- Transported label
        label& k()
        {
            return k_;
        }



        // Tracking

            //- Track all particles to their end point
            bool move(Cloud<trackedParticle>&, trackingData&, const scalar);

            //- Overridable function to handle the particle hitting a patch
            //  Executed before other patch-hitting functions
            bool hitPatch(Cloud<trackedParticle>&, trackingData&);

            //- Overridable function to handle the particle hitting a wedge
            void hitWedgePatch(Cloud<trackedParticle>&, trackingData&);

            //- Overridable function to handle the particle hitting a
            //  symmetry plane
            void hitSymmetryPlanePatch(Cloud<trackedParticle>&, trackingData&);

            //- Overridable function to handle the particle hitting a
            //  symmetry patch
            void hitSymmetryPatch(Cloud<trackedParticle>&, trackingData&);

            //- Overridable function to handle the particle hitting a cyclic
            void hitCyclicPatch(Cloud<trackedParticle>&, trackingData&);

            //- Overridable function to handle the particle hitting a cyclicAMI
            void hitCyclicAMIPatch
            (
                const vector&,
                const scalar,
                Cloud<trackedParticle>&,
                trackingData&
            );

            //- Overridable function to handle the particle hitting a cyclicACMI
            void hitCyclicACMIPatch
            (
                const vector&,
                const scalar,
                Cloud<trackedParticle>&,
                trackingData&
            );

            //- Overridable function to handle the particle hitting a
            //  cyclicRepeatAMI
            void hitCyclicRepeatAMIPatch
            (
                const vector&,
                const scalar,
                Cloud<trackedParticle>&,
                trackingData&
            );

            //- Overridable function to handle the particle hitting a
            //- processorPatch
            void hitProcessorPatch(Cloud<trackedParticle>&, trackingData&);

            //- Overridable function to handle the particle hitting a wallPatch
            void hitWallPatch(Cloud<trackedParticle>&, trackingData&);

            //- Convert processor patch addressing to the global equivalents
            //  and set the celli to the face-neighbour
            void correctAfterParallelTransfer(const label, trackingData&);


    // Ostream Operator

        friend Ostream& operator<<(Ostream&, const trackedParticle&);
};


template<>
inline bool contiguous<trackedParticle>()
{
    return true;
}

//template<>
//void Cloud<trackedParticle>::readFields();
//
//template<>
//void Cloud<trackedParticle>::writeFields() const;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
