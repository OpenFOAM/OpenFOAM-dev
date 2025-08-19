/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2025 OpenFOAM Foundation
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
    Foam::solidParticle

Description
    Simple solid spherical particle class with one-way coupling with the
    continuous phase.

SourceFiles
    solidParticleI.H
    solidParticle.C
    solidParticleIO.C

\*---------------------------------------------------------------------------*/

#ifndef solidParticle_H
#define solidParticle_H

#include "particle.H"
#include "Cloud.H"
#include "IOstream.H"
#include "autoPtr.H"
#include "cellPoint_interpolation.H"
#include "contiguous.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class solidParticleCloud;


// Forward declaration of friend functions and operators

class solidParticle;

Ostream& operator<<(Ostream&, const solidParticle&);


/*---------------------------------------------------------------------------*\
                        Class solidParticle Declaration
\*---------------------------------------------------------------------------*/

class solidParticle
:
    public particle
{
    // Private Data

        //- Size in bytes of the fields
        static const std::size_t sizeofFields_;

        //- Diameter
        scalar d_;

        //- Velocity of parcel
        vector U_;


public:

    friend class lagrangian::Cloud<solidParticle>;


    //- Class used to pass tracking data to the trackToFace function
    class trackingData
    :
        public particle::trackingData
    {
        // Interpolators for continuous phase fields

            const interpolations::cellPoint<scalar>& rhoInterp_;
            const interpolations::cellPoint<vector>& UInterp_;
            const interpolations::cellPoint<scalar>& nuInterp_;

        //- Local gravitational or other body-force acceleration
        const vector& g_;


    public:

        // Constructors

            inline trackingData
            (
                const solidParticleCloud& spc,
                const interpolations::cellPoint<scalar>& rhoInterp,
                const interpolations::cellPoint<vector>& UInterp,
                const interpolations::cellPoint<scalar>& nuInterp,
                const vector& g
            );


        // Member Functions

            inline const interpolations::cellPoint<scalar>& rhoInterp() const;

            inline const interpolations::cellPoint<vector>& UInterp() const;

            inline const interpolations::cellPoint<scalar>& nuInterp() const;

            inline const vector& g() const;
    };


    //- Runtime type information
    TypeName("solidParticle");


    // Constructors

        //- Construct from Istream
        solidParticle(Istream& is, bool readFields = true);

        //- Construct and return a clone
        virtual autoPtr<particle> clone() const
        {
            return autoPtr<particle>(new solidParticle(*this));
        }

        //- Construct from Istream and return
        static autoPtr<solidParticle> New(Istream& is)
        {
            return autoPtr<solidParticle>(new solidParticle(is));
        }


    // Member Functions

        // Access

            //- Return diameter
            inline scalar d() const;

            //- Return velocity
            inline const vector& U() const;


        // Tracking

            //- Move
            bool move(solidParticleCloud&, trackingData&);


        // Patch interactions

            //- Overridable function to handle the particle hitting a wallPatch
            void hitWallPatch(solidParticleCloud& cloud, trackingData& td);

            //- Transform the physical properties of the particle
            //  according to the given transformation
            virtual void transformProperties(const transformer&);


    // I-O

        static void readFields(lagrangian::Cloud<solidParticle>& c);

        static void writeFields(const lagrangian::Cloud<solidParticle>& c);


    // Ostream Operator

        friend Ostream& operator<<(Ostream&, const solidParticle&);
};


template<>
inline bool contiguous<solidParticle>()
{
    return true;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "solidParticleI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
