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
    Foam::IOList

Description
    A List of objects of type \<T\> with automated input and output.

SourceFiles
    IOList.C

\*---------------------------------------------------------------------------*/

#ifndef IOList_H
#define IOList_H

#include "List.H"
#include "regIOobject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class IOList Declaration
\*---------------------------------------------------------------------------*/

template<class T>
class IOList
:
    public regIOobject,
    public List<T>
{

public:

    //- Runtime type information
    TypeName("List");


    // Constructors

        //- Construct from IOobject
        IOList(const IOobject&);

        //- Construct from IOobject and size of IOList
        IOList(const IOobject&, const label);

        //- Construct from IOobject and a List
        IOList(const IOobject&, const List<T>&);

        //- Construct by transferring the List contents
        IOList(const IOobject&, const Xfer<List<T> >&);


    //- Destructor
    virtual ~IOList();


    // Member functions

        bool writeData(Ostream&) const;


    // Member operators

        void operator=(const IOList<T>&);

        void operator=(const List<T>&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "IOList.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
