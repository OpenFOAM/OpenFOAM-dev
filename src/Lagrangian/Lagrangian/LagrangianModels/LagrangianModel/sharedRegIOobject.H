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
    Foam::sharedRegIOobject

Description
    Container for a pointer to a shared object, who's memory is managed by an
    object registry. Object lifetime is managed by reference counting, which
    limits the object types to those that inherit a reference count.

SourceFiles
    sharedRegIOobjectI.H

\*---------------------------------------------------------------------------*/

#ifndef sharedRegIOobject_H
#define sharedRegIOobject_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class IOobject;

/*---------------------------------------------------------------------------*\
                          Class sharedRegIOobject Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class sharedRegIOobject
{
    // Private Data

        //- Pointer to the registered object
        Type* ptr_;


public:

    // Constructors

        //- Default construct
        inline sharedRegIOobject();

        //- Construct from an IO object and additional arguments
        template<class ... Args>
        inline sharedRegIOobject(const IOobject& io, const Args& ... args);

        //- Disallow default bitwise copy construction
        sharedRegIOobject(const sharedRegIOobject<Type>&) = delete;


    //- Destructor
    ~sharedRegIOobject();


    // Member Functions

        //- Is this object set?
        inline bool valid() const;

        //- Set the object
        template<class ... Args>
        inline void set(const IOobject& io, const Args& ... args);

        //- Return a reference to the object
        inline Type& ref() const;

        //- Clear the object
        inline void clear();


    // Member Operators

        //- Return a reference to the object
        inline const Type& operator()() const;

        //- Disallow default bitwise assignment
        void operator=(const sharedRegIOobject<Type>&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "sharedRegIOobjectI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
