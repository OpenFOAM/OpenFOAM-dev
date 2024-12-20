/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2024 OpenFOAM Foundation
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
    Foam::restartableRandomGenerator

Description
    Random number generator with the additional ability to go back to an
    earlier stored state. Useful for processes that occur within converging
    iteration loops to ensure that the random sequence is the same for each
    iteration, and is not preventing convergence by "re-randomising" the
    solution on each iteration.

\*---------------------------------------------------------------------------*/

#ifndef restartableRandomGenerator_H
#define restartableRandomGenerator_H

#include "randomGenerator.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class restartableRandomGenerator Declaration
\*---------------------------------------------------------------------------*/

class restartableRandomGenerator
:
    public randomGenerator
{
    // Private Data

        //- The random generator state at the start of the current sequence
        randomGenerator startRndGen_;


public:

    // Constructors

        //- Forward all constructors
        template<class ... Args>
        inline restartableRandomGenerator(Args&& ... args)
        :
            randomGenerator(std::forward<Args>(args) ...),
            startRndGen_(static_cast<const randomGenerator&>(*this))
        {}


    // Member Functions

        //- Start a sequence of random numbers. Either a new sequence, or a
        //  repeat of the previous, as indicated by the provided flag.
        inline void start(const bool repeat)
        {
            repeat
          ? static_cast<randomGenerator&>(*this) = startRndGen_
          : startRndGen_ = static_cast<randomGenerator&>(*this);
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#endif

// ************************************************************************* //
