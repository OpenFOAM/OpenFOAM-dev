/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2023 OpenFOAM Foundation
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
    Foam::LPtrList

Description
    Template class for non-intrusive linked PtrLists.

SourceFiles
    LPtrList.C
    LPtrListIO.C

\*---------------------------------------------------------------------------*/

#ifndef LPtrList_H
#define LPtrList_H

#include "ULPtrList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

template<class LListBase, class T> class LPtrList;

template<class LListBase, class T>
Istream& operator>>
(
    Istream&,
    LPtrList<LListBase, T>&
);


/*---------------------------------------------------------------------------*\
                          Class LPtrList Declaration
\*---------------------------------------------------------------------------*/

template<class LListBase, class T>
class LPtrList
:
    public ULPtrList<LListBase, T>
{
    // Private Member Functions

        //- Read from Istream using given Istream constructor class
        template<class INew>
        void read(Istream&, const INew&);


public:

    // Constructors

        //- Null construct
        LPtrList()
        {}

        //- Construct given initial T
        LPtrList(T* a)
        :
            ULPtrList<LListBase, T>(a)
        {}

        //- Construct from Istream using given Istream constructor class
        template<class INew>
        LPtrList(Istream&, const INew&);

        //- Construct from Istream using default Istream constructor class
        LPtrList(Istream&);

        //- Copy constructor
        LPtrList(const LPtrList&);

        //- Move constructor
        LPtrList(LPtrList&&);


    //- Destructor
    ~LPtrList();


    // Member Functions

        // Edit

            //- Remove the head element from the list and delete the pointer
            bool eraseHead();

            //- Clear the contents of the list
            void clear();

            //- Transfer the contents of the argument into this List
            //  and annul the argument list.
            void transfer(LPtrList<LListBase, T>&);


    // Member Operators

        //- Assignment operator
        void operator=(const LPtrList<LListBase, T>&);

        //- Move assignment operator
        void operator=(LPtrList<LListBase, T>&&);


    // IOstream Operators

        friend Istream& operator>> <LListBase, T>
        (
            Istream&,
            LPtrList<LListBase, T>&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "LPtrList.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
