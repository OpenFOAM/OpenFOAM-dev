/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2025 OpenFOAM Foundation
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
    Foam::oneOrTmp

Description
    Variant type that either contains and behaves as tmp object or represents a
    value of one. Provides multiplication and scaling operations which are
    optimised away if representing a value of one.

SourceFiles
    oneOrTmpI.H

\*---------------------------------------------------------------------------*/

#ifndef oneOrTmp_H
#define oneOrTmp_H

#include "one.H"
#include "tmp.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class oneOrTmp Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class oneOrTmp
{
    // Private Data

        //- Is this a value of one?
        bool isOne_;

        //- Temporary object
        tmp<Type> tt_;


public:

    // Constructors

        //- Construct as one
        inline oneOrTmp();

        //- Construct from an object
        inline oneOrTmp(const Type&);

        //- Construct from a temporary object
        inline oneOrTmp(const tmp<Type>&);

        //- Copy construct
        inline oneOrTmp(const oneOrTmp<Type>&);


    //- Destructor
    ~oneOrTmp();


    // Member Functions

        //- Scale an object by this value
        template<class Other>
        void scale(Other& o) const;


    // Member Operators

        //- Assign to one
        inline void operator=(const one&);

        //- Assign to a temporary object
        inline void operator=(const tmp<Type>&);

        //- Return an instance of the class which refers to this one. Used in
        //  the same way as tmp::operator() to prevent the tmp being consumed.
        inline oneOrTmp<Type> operator()() const;

        //- Multiply with an object
        template<class Other>
        inline tmp<Other> operator*(const Other&);

        //- Multiply with a temporary object
        template<class Other>
        inline tmp<Other> operator*(const tmp<Other>&);

        //- Disallow default bitwise assignment
        void operator=(const oneOrTmp<Type>&) = delete;
};


// Global Operators

//- Multiply with an object
template<class Type, class Other>
inline tmp<Other> operator*(const Other&, const oneOrTmp<Type>&);

//- Multiply with a temporary object
template<class Type, class Other>
inline tmp<Other> operator*(const tmp<Other>&, const oneOrTmp<Type>&);

//- Multiply-assign with an object
template<class Type, class Other>
inline void operator*=(Other&, const oneOrTmp<Type>&);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "oneOrTmpI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
