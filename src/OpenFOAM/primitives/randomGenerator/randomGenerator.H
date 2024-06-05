/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2018-2024 OpenFOAM Foundation
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
    Foam::randomGenerator

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
    randomGenerator.C
    randomGeneratorI.H

\*---------------------------------------------------------------------------*/

#ifndef randomGenerator_H
#define randomGenerator_H

#include "scalarField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class Istream;
class Ostream;

// Forward declaration of friend functions and operators
class randomGenerator;
Istream& operator>>(Istream&, randomGenerator&);
Ostream& operator<<(Ostream&, const randomGenerator&);

/*---------------------------------------------------------------------------*\
                       Class randomGenerator Declaration
\*---------------------------------------------------------------------------*/

class randomGenerator
{
public:

    // Public Classes

        //- Seed class
        class seed
        {
            // Private Data

                //- The seed value
                const uint64_t s_;


            // Private Member Functions

                //- Return the initial integer
                inline uint64_t x(const bool global) const;


        public:

            //- Allow randomGenerator to access the seed value
            friend class randomGenerator;


            // Constructors

                //- Construct from a label
                inline seed(const label s);

                //- Construct from a word
                inline seed(const word& s);
        };


private:

    // Private Static Data

        //- The parameters of the linear congruential iteration
        static const uint64_t A = 0x5DEECE66D, C = 0xB, M = uint64_t(1) << 48;


    // Private Data

        //- Is this generator global?
        const bool global_;

        //- The stored integer
        uint64_t x_;


    // Private Member Functions

        //- Check the state is synchronised
        void checkSync() const;

        //- Advance the state and return an integer sample
        inline uint64_t sample();


        // Scalars

            //- Return a scalar uniformly distributed between zero and one.
            //  Don't check synchronisation.
            inline scalar scalar01NoCheckSync();

            //- Return a scalar uniformly distributed between zero and one.
            //  Don't check synchronisation.
            inline scalar scalarABNoCheckSync(const scalar a, const scalar b);


        // Other Types

            //- Return a type with components uniformly distributed between
            //  zero and one. Don't check synchronisation.
            template<class Type>
            inline Type sample01NoCheckSync();

            //- Return a type with components uniformly distributed between two
            //  limits. Don't check synchronisation.
            template<class Type>
            inline Type sampleABNoCheckSync(const Type& a, const Type& b);


public:

    // Constructors

        //- Construct from a seed
        inline randomGenerator(const seed s, const bool global = false);

        //- Copy construct
        inline randomGenerator(const randomGenerator&);

        //- Move construct
        inline randomGenerator(randomGenerator&&);

        //- Construct from a stream
        randomGenerator(Istream& is, const bool global = false);

        //- Construct from a dictionary
        randomGenerator
        (
            const word& name,
            const dictionary& dict,
            randomGenerator&& defaultRndGen
        );

        //- Construct from a dictionary
        randomGenerator
        (
            const word& name,
            const dictionary& dict,
            const seed defaultS,
            const bool global = false
        );


    //- Destructor
    inline ~randomGenerator();


    // Member Functions

        // Scalars

            //- Return a scalar uniformly distributed between zero and one
            inline scalar scalar01();

            //- Return scalars uniformly distributed between zero and one
            inline tmp<scalarField> scalar01(const label n);

            //- Return a scalar uniformly distributed between two limits
            inline scalar scalarAB(const scalar a, const scalar b);

            //- Return scalars uniformly distributed between two limits
            inline tmp<scalarField> scalarAB
            (
                const label n,
                const scalar a,
                const scalar b
            );


        // Other Types

            //- Return a type with components uniformly distributed between
            //  zero and one
            template<class Type>
            inline Type sample01();

            //- Return types with components uniformly distributed between zero
            //  and one
            template<class Type>
            inline tmp<Field<Type>> sample01(const label n);

            //- Return a type with components uniformly distributed between two
            //  limits
            template<class Type>
            inline Type sampleAB(const Type& a, const Type& b);

            //- Return types with components uniformly distributed between two
            //  limits
            template<class Type>
            inline tmp<Field<Type>> sampleAB
            (
                const label n,
                const Type& a,
                const Type& b
            );


        //- Randomly permute a list
        template<class Container>
        inline void permute(Container& l);

        //- Create a randomly seeded generator
        inline randomGenerator generator();


    // Member Operators

        //- Copy-assignment
        void operator=(const randomGenerator&);

        //- Move-assignment
        void operator=(randomGenerator&&);


    // IOstream Operators

        //- Read from stream
        friend Istream& operator>>(Istream&, randomGenerator&);

        //- Write to stream
        friend Ostream& operator<<(Ostream&, const randomGenerator&);
};


template<>
inline scalar randomGenerator::sample01NoCheckSync();

template<>
inline label randomGenerator::sample01NoCheckSync();

template<>
inline scalar randomGenerator::sampleABNoCheckSync
(
    const scalar& a,
    const scalar& b
);

template<>
inline label randomGenerator::sampleABNoCheckSync
(
    const label& a,
    const label& b
);


void writeEntry(Ostream& os, const randomGenerator& rndGen);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "randomGeneratorI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
