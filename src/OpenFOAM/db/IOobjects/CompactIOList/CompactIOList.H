/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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

#include "ListCompactIO.H"
#include "regIOobject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class CompactIOList Declaration
\*---------------------------------------------------------------------------*/

template<class T, class BaseType>
class CompactIOList
:
    public regIOobject,
    public ListCompactIO<T, BaseType>
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

        //- Move construct by transferring the List contents
        CompactIOList(const IOobject&, List<T>&&);

        //- Move constructor
        CompactIOList(CompactIOList<T, BaseType>&&);


    // Destructor

        virtual ~CompactIOList();


    // Member Functions

        virtual bool writeObject
        (
            IOstream::streamFormat,
            IOstream::versionNumber,
            IOstream::compressionType,
            const bool write
        ) const;

        virtual bool writeData(Ostream&) const;


    // Member Operators

        void operator=(const CompactIOList<T, BaseType>&);
        void operator=(CompactIOList<T, BaseType>&&);

        void operator=(const List<T>&);
        void operator=(List<T>&&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "CompactIOList.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
