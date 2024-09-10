/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2024 OpenFOAM Foundation
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

#include "Function1.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class TableReader Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class TableReader
{
protected:

    // Protected Member Functions

        //- Do unit conversions on input data
        void convertRead
        (
            const Function1s::unitConversions& units,
            List<Tuple2<scalar, Type>>& table
        ) const;

        //- Do unit conversions on input data
        List<Tuple2<scalar, Type>> convertRead
        (
            const Function1s::unitConversions& units,
            const List<Tuple2<scalar, Type>>& table
        ) const;

        //- Do unit conversions in order to output data
        List<Tuple2<scalar, Type>> convertWrite
        (
            const Function1s::unitConversions& units,
            const List<Tuple2<scalar, Type>>& table
        ) const;


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
                const Function1s::unitConversions& units,
                const dictionary& dict
            ),
            (name, units, dict)
        );


    // Constructors

        //- Default construct
        TableReader();

        //- Construct and return a clone
        virtual autoPtr<TableReader<Type>> clone() const = 0;


    // Selector
    static autoPtr<TableReader<Type>> New
    (
        const word& name,
        const Function1s::unitConversions& units,
        const dictionary& dict
    );


    //- Destructor
    virtual ~TableReader();


    // Member Functions

        //- Read values
        virtual List<Tuple2<scalar, Type>> read
        (
            const Function1s::unitConversions& units,
            const dictionary& dict,
            const word& valuesKeyword="values"
        ) const = 0;

        //- Write settings and values
        virtual void write
        (
            Ostream& os,
            const Function1s::unitConversions& units,
            const List<Tuple2<scalar, Type>>& table,
            const word& valuesKeyword="values"
        ) const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define defineTableReader(Type)                                                \
                                                                               \
    defineNamedTemplateTypeNameAndDebug(TableReader<Type>, 0);                 \
                                                                               \
    defineTemplateRunTimeSelectionTable(TableReader<Type>, dictionary);


#define addTableReader(SS, Type)                                               \
                                                                               \
    defineNamedTemplateTypeNameAndDebug(TableReaders::SS<Type>, 0);            \
                                                                               \
    addTemplatedToRunTimeSelectionTable(TableReader, SS, Type, dictionary)


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "TableReader.C"
    #include "TableReaderNew.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
