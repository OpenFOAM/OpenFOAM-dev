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
    Foam::CompactIOList

Description
    A List of objects of type \<T\> with automated input and output using
    a compact storage. Behaves like IOList except when binary output in
    case it writes a CompactListList.

    Useful for lists of small sublists e.g. faceList, cellList.

SourceFiles
    CompactIOList.C

\*---------------------------------------------------------------------------*/

#ifndef CompactIOList_H
#define CompactIOList_H

#include "IOList.H"
#include "regIOobject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class Istream;
class Ostream;

// Forward declaration of friend functions and operators
template<class T, class BaseType> class CompactIOList;

template<class T, class BaseType> Istream& operator>>
(
    Istream&,
    CompactIOList<T, BaseType>&
);
template<class T, class BaseType> Ostream& operator<<
(
    Ostream&,
    const CompactIOList<T, BaseType>&
);

/*---------------------------------------------------------------------------*\
                           Class CompactIOList Declaration
\*---------------------------------------------------------------------------*/

template<class T, class BaseType>
class CompactIOList
:
    public regIOobject,
    public List<T>
{
    // Private Member Functions

        //- Read according to header type
        void readFromStream();

public:

    //- Runtime type information
    TypeName("CompactList");


    // Constructors

        //- Construct from IOobject
        CompactIOList(const IOobject&);

        //- Construct from IOobject and size of CompactIOList
        CompactIOList(const IOobject&, const label);

        //- Construct from IOobject and a List
        CompactIOList(const IOobject&, const List<T>&);

        //- Construct by transferring the List contents
        CompactIOList(const IOobject&, const Xfer<List<T> >&);


    // Destructor

        virtual ~CompactIOList();


    // Member functions

        virtual bool writeObject
        (
            IOstream::streamFormat,
            IOstream::versionNumber,
            IOstream::compressionType
        ) const;

        virtual bool writeData(Ostream&) const;


    // Member operators

        void operator=(const CompactIOList<T, BaseType>&);

        void operator=(const List<T>&);


    // IOstream operators

        //- Read List from Istream, discarding contents of existing List.
        friend Istream& operator>> <T, BaseType>
        (
            Istream&,
            CompactIOList<T, BaseType>&
        );

        // Write List to Ostream.
        friend Ostream& operator<< <T, BaseType>
        (
            Ostream&,
            const CompactIOList<T, BaseType>&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "CompactIOList.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
