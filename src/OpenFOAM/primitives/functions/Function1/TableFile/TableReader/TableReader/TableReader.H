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
    Foam::TableReader

Description
    Base class to read table data for tables

SourceFiles
    TableReader.C

\*---------------------------------------------------------------------------*/

#ifndef TableReader_H
#define TableReader_H

#include "dictionary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class TableReader Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class TableReader
{

public:

    //- Runtime type information
    TypeName("TableReader");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            TableReader,
            dictionary,
            (
                const word& name,
                const dictionary& dict,
                List<Tuple2<scalar, Type>>& table
            ),
            (name, dict, table)
        );


    // Constructors

        //- Construct from dictionary
        TableReader(const dictionary& dict);

        //- Construct and return a clone
        virtual autoPtr<TableReader<Type>> clone() const = 0;


    // Selector
    static autoPtr<TableReader<Type>> New
    (
        const word& name,
        const dictionary& dict,
        List<Tuple2<scalar, Type>>& table
    );


    //- Destructor
    virtual ~TableReader();


    // Member Functions

        //- Write additional information
        virtual void write
        (
            Ostream& os,
            const List<Tuple2<scalar, Type>>& table
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define defineTableReader(Type)                                                \
                                                                               \
    defineNamedTemplateTypeNameAndDebug(TableReader<Type>, 0);                 \
                                                                               \
    defineTemplateRunTimeSelectionTable(TableReader<Type>, dictionary);


#define makeTableReader(SS, Type)                                              \
                                                                               \
    defineNamedTemplateTypeNameAndDebug(TableReaders::SS<Type>, 0);            \
                                                                               \
    TableReader<Type>::adddictionaryConstructorToTable<TableReaders::SS<Type>> \
        add##SS##Type##ConstructorToTable_;

#define makeTableReaders(Type)                                                 \
    defineTableReader(Type);                                                   \
    makeTableReader(Embedded, Type);                                           \
    makeTableReader(Foam, Type);                                               \
    makeTableReader(Csv, Type)


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "TableReader.C"
    #include "TableReaderNew.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
