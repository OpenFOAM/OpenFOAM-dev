/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2018-2022 OpenFOAM Foundation
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
    Foam::Random

Description
    Random number generator

    This is a clone of the drand48 algorithm. This is significantly quicker
    than drand48, presumably due to the compiler inlining the sampling methods.
    It is also significantly quicker than the standard library linear
    congruential engine, as it does not use Schrage's algorithm to prevent
    overflow.

    See <http://pubs.opengroup.org/onlinepubs/007908775/xsh/drand48.html> for
    details of the seeding and iteration sequence.

SourceFiles
    RandomI.H

\*---------------------------------------------------------------------------*/

#ifndef Random_H
#define Random_H

#include "scalar.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class Random Declaration
\*---------------------------------------------------------------------------*/

class Random
{
    // Private Typedefs

        //- Working type
        typedef uint64_t type;


    // Private static data

        //- The parameters of the linear congruential iteration
        static const type A = 0x5DEECE66D, C = 0xB, M = type(1) << 48;


    // Private Data

        //- The stored integer
        type x_;

        //- Is a normal scalar sample stored?
        bool scalarNormalStored_;

        //- A stored normal scalar sample
        scalar scalarNormalValue_;


    // Private Member Functions

        //- Advance the state and return an integer sample
        inline type sample();


public:

    // Constructors

        //- Construct from a seed
        inline Random(const label s);


    //- Destructor
    inline ~Random();


    // Member Functions

        // Scalars

            //- Advance the state and return a scalar sample from a uniform
            //  distribution between zero and one
            inline scalar scalar01();

            //- Advance the state and return a scalar sample from a uniform
            //  distribution between two limits
            inline scalar scalarAB(const scalar a, const scalar b);

            //- Advance the state and return a scalar sample from a normal
            //  distribution with mean zero and standard deviation one
            scalar scalarNormal();


        // Other types

            //- Advance the state and return a sample of a given type from a
            //  uniform distribution between zero and one
            template<class Type>
            inline Type sample01();

            //- Advance the state and return a sample of a given type from a
            //  uniform distribution between two limits
            template<class Type>
            inline Type sampleAB(const Type& a, const Type& b);

            //- Advance the state and return a sample of a given type from a
            //  normal distribution with mean zero and standard deviation one
            template<class Type>
            inline Type sampleNormal();


        // Global scalars

            //- Advance the state and return a scalar sample from a uniform
            //  distribution between zero and one. Synchronises across all
            //  cores. Use of this is discouraged. It is expensive and
            //  introduces non-randomness in all cores other then the master.
            scalar globalScalar01();


        //- Randomly permute a list
        template<class Container>
        inline void permute(Container& l);
};


template<>
inline scalar Random::sample01();

template<>
inline label Random::sample01();

template<>
inline scalar Random::sampleAB(const scalar& a, const scalar& b);

template<>
inline label Random::sampleAB(const label& a, const label& b);

template<>
inline scalar Random::sampleNormal();


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "RandomI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
