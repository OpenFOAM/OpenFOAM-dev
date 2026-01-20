/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2026 OpenFOAM Foundation
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

template<class Coordinate, class Value>
class TableReader
{
protected:

    // Protected Member Functions

        //- Do unit conversions on input data
        void convertRead
        (
            const Function1s::unitConversions& units,
            List<Tuple2<Coordinate, Value>>& table
        ) const;

        //- Do unit conversions on input data
        List<Tuple2<Coordinate, Value>> convertRead
        (
            const Function1s::unitConversions& units,
            const List<Tuple2<Coordinate, Value>>& table
        ) const;

        //- Do unit conversions in order to output data
        List<Tuple2<Coordinate, Value>> convertWrite
        (
            const Function1s::unitConversions& units,
            const List<Tuple2<Coordinate, Value>>& table
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
        virtual autoPtr<TableReader<Coordinate, Value>> clone() const = 0;


    // Selector
    static autoPtr<TableReader<Coordinate, Value>> New
    (
        const word& name,
        const Function1s::unitConversions& units,
        const dictionary& dict
    );


    //- Destructor
    virtual ~TableReader();


    // Member Functions

        //- Read values
        virtual List<Tuple2<Coordinate, Value>> read
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
            const List<Tuple2<Coordinate, Value>>& table,
            const word& valuesKeyword="values"
        ) const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define defineTableReader(Coordinate, Value)                                   \
                                                                               \
    typedef TableReader<Coordinate, Value>                                     \
        TableReader##Coordinate##Value;                                        \
                                                                               \
    defineNamedTemplateTypeNameAndDebug                                        \
    (                                                                          \
        TableReader##Coordinate##Value,                                        \
        0                                                                      \
    );                                                                         \
                                                                               \
    defineTemplateRunTimeSelectionTable                                        \
    (                                                                          \
        TableReader##Coordinate##Value,                                        \
        dictionary                                                             \
    );


#define addTableReader(TableReaderType, Coordinate, Value)                     \
                                                                               \
    typedef TableReader<Coordinate, Value>                                     \
        TableReader##Coordinate##Value;                                        \
                                                                               \
    namespace TableReaders                                                     \
    {                                                                          \
        typedef TableReaderType<Coordinate, Value>                             \
            TableReaderType##TableReader##Coordinate##Value;                   \
                                                                               \
        defineNamedTemplateTypeNameAndDebug                                    \
        (                                                                      \
            TableReaderType##TableReader##Coordinate##Value,                   \
            0                                                                  \
        );                                                                     \
                                                                               \
        addToRunTimeSelectionTable                                             \
        (                                                                      \
            TableReader##Coordinate##Value,                                    \
            TableReaderType##TableReader##Coordinate##Value,                   \
            dictionary                                                         \
        );                                                                     \
    }


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "TableReader.C"
    #include "TableReaderNew.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
