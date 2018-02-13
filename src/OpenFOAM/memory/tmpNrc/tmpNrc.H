/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016-2018 OpenFOAM Foundation
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
    Foam::tmpNrc

Description
    A class for managing temporary objects without reference counting.

SourceFiles
    tmpNrcI.H

See also
    Foam::autoPtr

\*---------------------------------------------------------------------------*/

#ifndef tmpNrc_H
#define tmpNrc_H

#include "word.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                             Class tmpNrc Declaration
\*---------------------------------------------------------------------------*/

template<class T>
class tmpNrc
{
    // Private data

        //- Object types
        enum type
        {
            TMP,
            CONST_REF
        };

        //- Type of object
        mutable type type_;

        //- Pointer to object
        mutable T* ptr_;


public:

    // Null reference counter class
    class refCount
    {
    public:

        refCount()
        {}
    };


    // Constructors

        //- Store object pointer
        inline explicit tmpNrc(T* = 0);

        //- Store object const reference
        inline tmpNrc(const T&);

        //- Construct copy and increment reference count
        inline tmpNrc(const tmpNrc<T>&);

        //- Construct copy transferring content of temporary if required
        inline tmpNrc(const tmpNrc<T>&, bool allowTransfer);


    //- Destructor: deletes temporary object when the reference count is 0
    inline ~tmpNrc();


    // Member Functions

        // Access

            //- Return true if this is really a temporary object
            inline bool isTmp() const;

            //- Return true if this temporary object empty,
            //  ie, a temporary without allocation
            inline bool empty() const;

            //- Is this temporary object valid,
            //  ie, it is a reference or a temporary that has been allocated
            inline bool valid() const;

            //- Return the type name of the tmpNrc
            //  constructed from the type name of T
            inline word typeName() const;


        // Edit

            //- Return non-const reference or generate a fatal error
            //  if the object is const.
            inline T& ref();

            //- Return tmpNrc pointer for reuse.
            //  Returns a clone if the object is not a temporary
            inline T* ptr() const;

            //- If object pointer points to valid object:
            //  delete object and set pointer to nullptr
            inline void clear() const;


    // Member operators

        //- Const dereference operator
        inline const T& operator()() const;

        //- Const cast to the underlying type reference
        inline operator const T&() const;

        //- Return object pointer
        inline T* operator->();

        //- Return const object pointer
        inline const T* operator->() const;

        //- Assignment to pointer changing this tmpNrc to a temporary T
        inline void operator=(T*);

        //- Assignment transferring the temporary T to this tmpNrc
        inline void operator=(const tmpNrc<T>&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "tmpNrcI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
