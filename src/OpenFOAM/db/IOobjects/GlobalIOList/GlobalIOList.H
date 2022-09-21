/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2015-2022 OpenFOAM Foundation
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
    Foam::GlobalIOList

Description
    IOContainer with global data (so optionally read from master)

SourceFiles
    GlobalIOList.C

\*---------------------------------------------------------------------------*/

#ifndef GlobalIOList_H
#define GlobalIOList_H

#include "regIOobject.H"
#include "List.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class GlobalIOListBase Declaration
\*---------------------------------------------------------------------------*/

template
<
    template<class> class Container,
    template<class> class IOContainer,
    class Type
>
class GlobalIOListBase
:
    public regIOobject,
    public Container<Type>
{

public:

    // Constructors

        //- Construct from IOobject
        GlobalIOListBase(const IOobject&);

        //- Construct from IOobject and size (does not set values)
        GlobalIOListBase(const IOobject&, const label size);

        //- Construct from IOobject and a Container
        GlobalIOListBase(const IOobject&, const Container<Type>&);

        //- Move construct by transferring the Container contents
        GlobalIOListBase(const IOobject&, Container<Type>&&);

        //- Copy constructor
        GlobalIOListBase(const GlobalIOListBase<Container, IOContainer, Type>&);

        //- Move constructor
        GlobalIOListBase(GlobalIOListBase<Container, IOContainer, Type>&&);


    //- Destructor
    virtual ~GlobalIOListBase();


    // Member Functions

        //- Is object global
        virtual bool global() const
        {
            return true;
        }

        //- ReadData function required for regIOobject read operation
        virtual bool readData(Istream&);

        //- WriteData function required for regIOobject write operation
        virtual bool writeData(Ostream&) const;

    // Member Operators

        //- Inherit assignment operators
        using Container<Type>::operator=;

        //- Assign to another GlobalIOList
        void operator=(const GlobalIOListBase<Container, IOContainer, Type>&);

        //- Move-assign to another GlobalIOList
        void operator=(GlobalIOListBase<Container, IOContainer, Type>&&);
};


/*---------------------------------------------------------------------------*\
                           Class GlobalIOList Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class GlobalIOList
:
    public GlobalIOListBase<List, GlobalIOList, Type>
{

public:

    //- Runtime type information
    TypeName("List");


    // Constructors

        //- Inherit constructors
        using GlobalIOListBase<List, GlobalIOList, Type>::GlobalIOListBase;


    // Member Operators

        //- Inherit assignment operators
        using GlobalIOListBase<List, GlobalIOList, Type>::operator=;
};


//- Trait for obtaining global status
template<class Type>
struct typeGlobal<GlobalIOList<Type>>
{
    static const bool global = true;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "GlobalIOList.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
