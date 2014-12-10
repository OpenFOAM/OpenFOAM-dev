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
    Foam::cachedRandom

Description
    Random number generator.

    Pre-computes and caches samples on construction, so that when sample01()
    is called, the function simply returns the next (pre-computed) sample. On
    reaching the last sample, the sample sequence is repeated.

    Constructed using a seed and sample count. If the supplied count is
    negative, no caching is performed, and a new sample is generated on each
    call to sample01().

    Note: the copy constructor cannot be used if count = -1.

SourceFiles
    cachedRandomI.H
    cachedRandom.C
    cachedRandomTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef cachedRandom_H
#define cachedRandom_H

#include "scalarList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class cachedRandom;

/*---------------------------------------------------------------------------*\
                       Class cachedRandom Declaration
\*---------------------------------------------------------------------------*/

class cachedRandom
{
private:

    // Private data

        //- Initial random number seed
        label seed_;

        //- List of scalar samples
        scalarList samples_;

        //- Current sample marker
        label sampleI_;


    // Private Member Functions

        //- Returns the current sample
        scalar scalar01();


public:


    // Constructors

        //- Construct given seed and sample count
        cachedRandom(const label seed, const label count);

        //- Copy constructor with optional reset of sampleI
        cachedRandom(const cachedRandom& cr, const bool reset = false);


    // Destructor
    ~cachedRandom();


    // Member functions

        // Access

            //- Return const access to the initial random number seed
            inline label seed() const;

            //- Return const access to the list of samples
            inline const scalarList& samples() const;

            //- Return the current sample marker
            inline label sampleI() const;


        // Manipulation

            //- Return non-const access to the sample marker
            inline label& sampleI();


        // Evaluation

            // Random numbers
                //- Return a sample whose components lie in the range 0-1
                template<class Type>
                Type sample01();

                //- Return a sample between start and end
                template<class Type>
                Type position(const Type& start, const Type& end);

                //- Randomise value in the range 0-1
                template<class Type>
                void randomise01(Type& value);


            // Global random numbers - consistent across all processors

                //- Return a sample whose components lie in the range 0-1
                template<class Type>
                Type globalSample01();

                //- Return a sample between start and end
                template<class Type>
                Type globalPosition(const Type& start, const Type& end);

                //- Randomise value in the range 0-1
                template<class Type>
                void globalRandomise01(Type& value);


        // Operators

            //- Assignment operator
            void operator=(const cachedRandom& cr);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Template specialisations

template<>
label cachedRandom::sample01<label>();

template<>
scalar cachedRandom::sample01<scalar>();

template<>
label cachedRandom::position<label>(const label& start, const label& end);

template<>
scalar cachedRandom::position<scalar>
(
    const scalar& start,
    const scalar& end
);

template<>
label cachedRandom::globalSample01<label>();

template<>
scalar cachedRandom::globalSample01<scalar>();

template<>
label cachedRandom::globalPosition<label>(const label& start, const label& end);

template<>
scalar cachedRandom::globalPosition<scalar>
(
    const scalar& start,
    const scalar& end
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "cachedRandomI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "cachedRandomTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
