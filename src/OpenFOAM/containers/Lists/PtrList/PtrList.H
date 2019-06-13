/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
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
    Foam::PtrList

Description
    A templated 1D list of pointers to objects of type \<T\>, where the
    size of the array is known and used for subscript bounds checking, etc.

    The element operator [] returns a reference to the object rather than a
    pointer.

SourceFiles
    PtrListI.H
    PtrList.C
    PtrListIO.C

\*---------------------------------------------------------------------------*/

#ifndef PtrList_H
#define PtrList_H

#include "UPtrList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes

template<class T> class autoPtr;
template<class T> class tmp;

class SLListBase;
template<class LListBase, class T> class LPtrList;
template<class T>
using SLPtrList = LPtrList<SLListBase, T>;

// Forward declaration of friend functions and operators
template<class T> class PtrList;

template<class T>
Istream& operator>>(Istream&, PtrList<T>&);


/*---------------------------------------------------------------------------*\
                           Class PtrList Declaration
\*---------------------------------------------------------------------------*/

template<class T>
class PtrList
:
    public UPtrList<T>
{

protected:

    // Protected Member Functions

        //- Read from Istream using given Istream constructor class
        template<class INew>
        void read(Istream&, const INew& inewt);


public:

    // Constructors

        //- Null Constructor
        PtrList();

        //- Construct with size specified
        explicit PtrList(const label);

        //- Copy constructor
        PtrList(const PtrList<T>&);

        //- Copy constructor with additional argument for clone
        template<class CloneArg>
        PtrList(const PtrList<T>&, const CloneArg&);

        //- Move constructor
        PtrList(PtrList<T>&&);

        //- Construct as copy or re-use as specified
        PtrList(PtrList<T>&, bool reuse);

        //- Construct as copy of SLPtrList<T>
        explicit PtrList(const SLPtrList<T>&);

        //- Construct from Istream using given Istream constructor class
        template<class INew>
        PtrList(Istream&, const INew&);

        //- Construct from Istream using default Istream constructor class
        PtrList(Istream&);


    //- Destructor
    ~PtrList();


    // Member Functions

        // Edit

            //- Reset size of PtrList. If extending the PtrList, new entries are
            //  set to nullptr. If truncating the PtrList, removed entries are
            //  deleted
            void setSize(const label);

            //- Alias for setSize(const label)
            inline void resize(const label);

            //- Clear the PtrList, i.e. set size to zero deleting all the
            //  allocated entries
            void clear();

            //- Append an element at the end of the list
            inline void append(T*);
            inline void append(const autoPtr<T>&);
            inline void append(const tmp<T>&);

            //- Transfer the contents of the argument PtrList into this PtrList
            //  and annul the argument list
            void transfer(PtrList<T>&);

            //- Is element set
            inline bool set(const label) const;

            //- Set element to given T* and return old element (can be nullptr)
            inline autoPtr<T> set(const label, T*);

            //- Set element to given autoPtr<T> and return old element
            inline autoPtr<T> set(const label, const autoPtr<T>&);

            //- Set element to given tmp<T> and return old element
            inline autoPtr<T> set(const label, const tmp<T>&);

            //- Reorders elements. Ordering does not have to be done in
            //  ascending or descending order. Reordering has to be unique.
            //  (is shuffle)
            void reorder(const labelUList& oldToNew);

            //- Reorders elements. Ordering does not have to be done in
            //  ascending or descending order. Reordering has to be unique.
            //  Note: can create unset elements
            void shuffle(const labelUList& newToOld);


    // Member Operators

        //- Assignment operator
        void operator=(const PtrList<T>&);

        //- Move assignment operator
        void operator=(PtrList<T>&&);


    // IOstream operator

        //- Read PtrList from Istream, discarding contents of existing PtrList
        friend Istream& operator>> <T>(Istream&, PtrList<T>&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "PtrListI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "PtrList.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
