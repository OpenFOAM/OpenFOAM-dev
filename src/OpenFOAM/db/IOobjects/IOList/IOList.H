/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2022 OpenFOAM Foundation
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
    A List of objects of type \<Type\> with automated input and output.

SourceFiles
    IOList.C

\*---------------------------------------------------------------------------*/

#ifndef IOList_H
#define IOList_H

#include "regIOobject.H"
#include "List.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class IOListBase Declaration
\*---------------------------------------------------------------------------*/

template
<
    template<class> class Container,
    template<class> class IOContainer,
    class Type
>
class IOListBase
:
    public regIOobject,
    public Container<Type>
{

public:

    // Public Static Data

        //- Does this object support re-reading?
        static const bool rereading = false;


    // Constructors

        //- Construct from IOobject
        IOListBase(const IOobject&);

        //- Construct from IOobject; does local processor require reading?
        IOListBase(const IOobject&, const bool read);

        //- Construct from IOobject and size (does not set values)
        IOListBase(const IOobject&, const label size);

        //- Construct from IOobject and a Container
        IOListBase(const IOobject&, const Container<Type>&);

        //- Move construct by transferring the Container contents
        IOListBase(const IOobject&, Container<Type>&&);

        //- Copy constructor
        IOListBase(const IOListBase<Container, IOContainer, Type>&);

        //- Move constructor
        IOListBase(IOListBase<Container, IOContainer, Type>&&);


    //- Destructor
    virtual ~IOListBase();


    // Member Functions

        //- WriteData function required for regIOobject write operation
        virtual bool writeData(Ostream&) const;


    // Member Operators

        //- Inherit assignment operators
        using Container<Type>::operator=;

        //- Assign to another IOList
        void operator=(const IOListBase<Container, IOContainer, Type>&);

        //- Move-assign to another IOList
        void operator=(IOListBase<Container, IOContainer, Type>&&);
};


/*---------------------------------------------------------------------------*\
                           Class IOList Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class IOList
:
    public IOListBase<List, IOList, Type>
{

public:

    //- Runtime type information
    TypeName("List");


    // Constructors

        //- Inherit constructors
        using IOListBase<List, IOList, Type>::IOListBase;


    // Member Operators

        //- Inherit assignment operators
        using IOListBase<List, IOList, Type>::operator=;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "IOList.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
