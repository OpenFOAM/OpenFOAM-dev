/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2013-2020 OpenFOAM Foundation
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
    Foam::findCellParticle

Description
    Particle class that finds cells by tracking

SourceFiles
    findCellParticle.C

\*---------------------------------------------------------------------------*/

#ifndef findCellParticle_H
#define findCellParticle_H

#include "particle.H"
#include "Cloud.H"
#include "autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class findCellParticleCloud;


// Forward declaration of friend functions and operators

class findCellParticle;

Ostream& operator<<(Ostream&, const findCellParticle&);


/*---------------------------------------------------------------------------*\
                      Class findCellParticle Declaration
\*---------------------------------------------------------------------------*/

class findCellParticle
:
    public particle
{
    // Private Data

        //- Displacement over which to track
        vector displacement_;

        //- Passive data
        label data_;


public:

    friend class Cloud<findCellParticle>;

    //- Class used to pass tracking data to the trackToFace function
    class trackingData
    :
        public particle::trackingData
    {
        labelListList& cellToData_;
        List<List<point>>& cellToEnd_;

    public:

        // Constructors

            trackingData
            (
                Cloud<findCellParticle>& cloud,
                labelListList& cellToData,
                List<List<point>>& cellToEnd
            )
            :
                particle::trackingData(cloud),
                cellToData_(cellToData),
                cellToEnd_(cellToEnd)
            {}


        // Member Functions

            labelListList& cellToData()
            {
                return cellToData_;
            }

            List<List<point>>& cellToEnd()
            {
                return cellToEnd_;
            }
    };


    // Constructors

        //- Construct from components
        findCellParticle
        (
            const polyMesh& mesh,
            const barycentric& coordinates,
            const label celli,
            const label tetFacei,
            const label tetPtI,
            const vector& displacement,
            const label data
        );

        //- Construct from a position and a cell, searching for the rest of the
        //  required topology
        findCellParticle
        (
            const polyMesh& mesh,
            const vector& position,
            const label celli,
            const vector& displacement,
            const label data
        );

        //- Construct from Istream
        findCellParticle
        (
            const polyMesh& mesh,
            Istream& is,
            bool readFields = true
        );

        //- Construct and return a clone
        autoPtr<particle> clone() const
        {
            return autoPtr<particle>(new findCellParticle(*this));
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

            autoPtr<findCellParticle> operator()(Istream& is) const
            {
                return autoPtr<findCellParticle>
                (
                    new findCellParticle(mesh_, is, true)
                );
            }
        };


    // Member Functions

        //- Displacement over which to track
        const vector& displacement() const
        {
            return displacement_;
        }

        //- Displacement over which to track
        vector& displacement()
        {
            return displacement_;
        }

        //- Transported label
        label data() const
        {
            return data_;
        }

        //- Transported label
        label& data()
        {
            return data_;
        }


        // Tracking

            //- Track all particles to their end point
            bool move(Cloud<findCellParticle>&, trackingData&, const scalar);

            //- Overridable function to handle the particle hitting a patch
            //  Executed before other patch-hitting functions
            bool hitPatch(Cloud<findCellParticle>&, trackingData&);

            //- Overridable function to handle the particle hitting a wedge
            void hitWedgePatch(Cloud<findCellParticle>&, trackingData&);

            //- Overridable function to handle the particle hitting a
            //  symmetry plane
            void hitSymmetryPlanePatch(Cloud<findCellParticle>&, trackingData&);

            //- Overridable function to handle the particle hitting a
            //  symmetry patch
            void hitSymmetryPatch(Cloud<findCellParticle>&, trackingData&);

            //- Overridable function to handle the particle hitting a cyclic
            void hitCyclicPatch(Cloud<findCellParticle>&, trackingData&);

            //- Overridable function to handle the particle hitting a cyclicAMI
            void hitCyclicAMIPatch
            (
                const vector&,
                const scalar,
                Cloud<findCellParticle>&,
                trackingData&
            );

            //- Overridable function to handle the particle hitting a cyclicACMI
            void hitCyclicACMIPatch
            (
                const vector&,
                const scalar,
                Cloud<findCellParticle>&,
                trackingData&
            );

            //- Overridable function to handle the particle hitting a
            //cyclicRepeatAMI
            void hitCyclicRepeatAMIPatch
            (
                const vector&,
                const scalar,
                Cloud<findCellParticle>&,
                trackingData&
            );

            //- Overridable function to handle the particle hitting a
            //- processorPatch
            void hitProcessorPatch(Cloud<findCellParticle>&, trackingData&);

            //- Overridable function to handle the particle hitting a wallPatch
            void hitWallPatch(Cloud<findCellParticle>&, trackingData&);


    // Ostream Operator

        friend Ostream& operator<<(Ostream&, const findCellParticle&);
};


template<>
inline bool contiguous<findCellParticle>()
{
    return true;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
