/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2017-2020 OpenFOAM Foundation
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
    Foam::Roots

Description
    Templated storage for the roots of polynomial equations, plus flags to
    indicate the nature of the roots.

SourceFiles
    RootsI.H
    Roots.C

\*---------------------------------------------------------------------------*/

#ifndef Roots_H
#define Roots_H

#include "VectorSpace.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- Types of root
enum class rootType
{
    real = 0,
    complex,
    posInf,
    negInf,
    nan
};


/*---------------------------------------------------------------------------*\
                            Class Roots Declaration
\*---------------------------------------------------------------------------*/

template<direction N>
class Roots
:
    public VectorSpace<Roots<N>, scalar, N>
{
    // Private Data

        //- Root types, encoded into a single integer
        label types_;

public:

    // Constructors

        //- Construct null
        inline Roots();

        //- Construct with a uniform value
        inline Roots(const rootType t, const scalar x);

        //- Construct by concatenation
        inline Roots
        (
            const rootType t,
            const scalar x,
            const Roots<N - 1>& xs
        );

        //- Construct by concatenation
        inline Roots
        (
            const Roots<N - 1>& xs,
            const rootType t,
            const scalar x
        );

        //- Construct by concatenation
        template <direction M>
        inline Roots(const Roots<M>& xs, const Roots<N - M>& ys);


    // Member Functions

        //- Set the type of the i-th root
        inline void type(const direction i, const rootType t);

        //- Return the type of the i-th root
        inline rootType type(const direction i) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "RootsI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
