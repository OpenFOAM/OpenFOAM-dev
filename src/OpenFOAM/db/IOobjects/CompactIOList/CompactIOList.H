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
    Foam::CompactIOList

Description
    A List of objects of type \<Type\> with automated input and output using
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

// Forward declaration of friend functions and operators
template
<
    template<class> class Container,
    template<class> class IOContainer,
    template<class> class CompactIOContainer,
    class Type
>
class CompactIOListBase;

template
<
    template<class> class Container,
    template<class> class IOContainer,
    template<class> class CompactIOContainer,
    class Type
>
void writeEntry
(
    Ostream&,
    const CompactIOListBase<Container, IOContainer, CompactIOContainer, Type>&
);

template
<
    template<class> class Container,
    template<class> class IOContainer,
    template<class> class CompactIOContainer,
    class Type
>
Istream& operator>>
(
    Istream&,
    CompactIOListBase<Container, IOContainer, CompactIOContainer, Type>&
);

template
<
    template<class> class Container,
    template<class> class IOContainer,
    template<class> class CompactIOContainer,
    class Type
>
Ostream& operator<<
(
    Ostream&,
    const CompactIOListBase<Container, IOContainer, CompactIOContainer, Type>&
);


/*---------------------------------------------------------------------------*\
                        Class CompactIOListBase Declaration
\*---------------------------------------------------------------------------*/

template
<
    template<class> class Container,
    template<class> class IOContainer,
    template<class> class CompactIOContainer,
    class Type
>
class CompactIOListBase
:
    public regIOobject,
    public Container<Type>
{
    // Private Member Functions

        //- Does the container have too many elements to represent with the
        //  current integer size type?
        bool overflows() const;

        //- Convert the data to a compact representation
        void convertToCompact
        (
            labelList& start,
            List<typename Type::value_type>& elems
        ) const;

        //- Convert the data from a compact representation
        void convertFromCompact
        (
            const labelList& start,
            const List<typename Type::value_type>& elems
        );

        //- Read according to header type
        void readFromStream(const bool read = true);


public:

    // Constructors

        //- Construct from IOobject
        CompactIOListBase(const IOobject&);

        //- Construct from IOobject; does local processor require reading?
        CompactIOListBase(const IOobject&, const bool read);

        //- Construct from IOobject and size of CompactIOListBase
        CompactIOListBase(const IOobject&, const label);

        //- Construct from IOobject and a Container
        CompactIOListBase(const IOobject&, const Container<Type>&);

        //- Move construct by transferring the Container contents
        CompactIOListBase(const IOobject&, Container<Type>&&);

        //- Move constructor
        CompactIOListBase
        (
            CompactIOListBase
            <
                Container,
                IOContainer,
                CompactIOContainer,
                Type
            >&&
        );


    //- Destructor
    virtual ~CompactIOListBase();


    // Member Functions

        //- WriteObject function. Changes the write type depending on whether
        //  the output is compact or not.
        virtual bool writeObject
        (
            IOstream::streamFormat,
            IOstream::versionNumber,
            IOstream::compressionType,
            const bool write
        ) const;

        //- WriteData function required for regIOobject write operation
        virtual bool writeData(Ostream&) const;


    // Member Operators

        //- Inherit assignment operators
        using Container<Type>::operator=;

         //- Assign to another CompactIOList
        void operator=
        (
            const CompactIOListBase
            <
                Container,
                IOContainer,
                CompactIOContainer,
                Type
            >&
        );

        //- Move-assign to another CompactIOList
        void operator=
        (
            CompactIOListBase
            <
                Container,
                IOContainer,
                CompactIOContainer,
                Type
            >&&
        );


    // IOstream functions

        //- Write list to Ostream
        friend void writeEntry
        <
            Container,
            IOContainer,
            CompactIOContainer,
            Type
        >
        (
            Ostream&,
            const CompactIOListBase
            <
                Container,
                IOContainer,
                CompactIOContainer,
                Type
            >&
        );


    // IOstream Operators

        //- Read list from Istream
        friend Istream& operator>>
        <
            Container,
            IOContainer,
            CompactIOContainer,
            Type
        >
        (
            Istream&,
            CompactIOListBase
            <
                Container,
                IOContainer,
                CompactIOContainer,
                Type
            >&
        );

        // Write list to Ostream
        friend Ostream& operator<<
        <
            Container,
            IOContainer,
            CompactIOContainer,
            Type
         >
        (
            Ostream&,
            const CompactIOListBase
            <
                Container,
                IOContainer,
                CompactIOContainer,
                Type
            >&
        );
};


/*---------------------------------------------------------------------------*\
                        Class CompactIOList Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class CompactIOList
:
    public CompactIOListBase<List, IOList, CompactIOList, Type>
{

public:

    //- Runtime type information
    TypeName("ListList");


    // Constructors

        //- Inherit constructors
        using
            CompactIOListBase<List, IOList, CompactIOList, Type>::
            CompactIOListBase;


    // Member Operators

        //- Inherit assignment operators
        using
            CompactIOListBase<List, IOList, CompactIOList, Type>::
            operator=;
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
