/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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
    Foam::Xfer

Description
    A simple container for copying or transferring objects of type \<T\>.

    The wrapped object of type \<T\> must implement a transfer() method and
    an operator=() copy method.

    Since it is decided upon construction of the Xfer object whether the
    parameter is to be copied or transferred, the contents of the resulting
    Xfer object can be transferred unconditionally. This greatly simplifies
    defining constructors or methods in other classes with mixed
    transfer/copy semantics without requiring 2^N different versions.

    When transferring between dissimilar types, the xferCopyTo() and
    xferMoveTo() functions can prove useful. An example is transferring
    from a DynamicList to a List.  Since the
    List\<T\>::transfer(List\<T\>&) method could result in some allocated
    memory becoming inaccessible, the xferMoveTo() function should be used to
    invoke the correct List\<T\>::transfer(DynamicList\<T\>&) method.

    \code
        DynamicList<label> dynLst;
        ...
        labelList plainLst( xferMoveTo<labelList>(dynLst) );
    \endcode

    Of course, since this example is a very common operation, the
    DynamicList::xfer() method transfers to a plain List anyhow.
    It would thus be simpler (and clearer) just to use the following code:

    \code
        DynamicList<label> dynLst;
        ...
        labelList plainLst(dynLst.xfer());
    \endcode

SeeAlso
    xferCopy, xferCopyTo, xferMove, xferMoveTo, xferTmp

SourceFiles
    XferI.H

\*---------------------------------------------------------------------------*/

#ifndef Xfer_H
#define Xfer_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
template<class T> class tmp;

/*---------------------------------------------------------------------------*\
                           Class Xfer Declaration
\*---------------------------------------------------------------------------*/

template<class T>
class Xfer
{
    // Private data

        //- Pointer to underlying datatype
        mutable T* ptr_;

public:

    // Constructors

        //- Store object pointer and manage its deletion
        //  Can also be used later to transfer by assignment
        inline explicit Xfer(T* = 0);

        //- Construct by copying or by transferring the parameter contents
        inline explicit Xfer(T&, bool allowTransfer=false);

        //- Construct by copying the parameter contents
        inline explicit Xfer(const T&);

        //- Construct by transferring the contents
        inline Xfer(const Xfer<T>&);

    //- Destructor
    inline ~Xfer();

    // Member Functions

        //- Return a null object reference
        inline static const Xfer<T>& null();

    // Member Operators

        //- Transfer the contents into the object
        inline void operator=(T&);

        //- Transfer the contents into the object
        inline void operator=(const Xfer<T>&);

        //- Reference to the underlying datatype
        inline T& operator()() const;

        //- Pointer to the underlying datatype
        inline T* operator->() const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

/**
 * Construct by copying the contents of the \a arg
 *
 * \sa xferCopyTo, xferMove, xferMoveTo, xferTmp and Foam::Xfer
*/
template<class T>
inline Xfer<T> xferCopy(const T&);

/**
 * Construct by transferring the contents of the \a arg
 *
 * \sa xferCopy, xferCopyTo, xferMoveTo, xferTmp and Foam::Xfer
*/
template<class T>
inline Xfer<T> xferMove(T&);


/**
 * Construct by transferring the contents of the \a arg
 *
 * \sa xferCopy, xferCopyTo, xferMove, xferMoveTo and Foam::Xfer
*/
template<class T>
inline Xfer<T> xferTmp(Foam::tmp<T>&);


/**
 * Construct by copying the contents of the \a arg
 * between dissimilar types
 *
 * \sa xferCopy, xferMove, xferMoveTo, xferTmp and Foam::Xfer
*/
template<class To, class From>
inline Xfer<To> xferCopyTo(const From&);


/**
 * Construct by transferring the contents of the \a arg
 * between dissimilar types
 *
 * \par Example Use
 * \code
 *     DynamicList<label> dynLst;
 *     ...
 *     labelList plainLst( xferMoveTo<labelList>(dynLst) );
 * \endcode
 *
 * \sa xferCopy, xferCopyTo, xferMove, xferTmp and Foam::Xfer
*/
template<class To, class From>
inline Xfer<To> xferMoveTo(From&);


} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "XferI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
