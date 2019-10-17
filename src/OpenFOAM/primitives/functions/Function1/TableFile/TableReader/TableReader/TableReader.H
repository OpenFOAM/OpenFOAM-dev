/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
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
#include "fieldTypes.H"
#include "runTimeSelectionTables.H"
#include "Tuple2.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class TableReader Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class TableReader
{
private:

    // Private Member Functions

        //- Read a 1D table
        virtual void read
        (
            ISstream&,
            List<Tuple2<scalar, Type>>&
        ) const = 0;

        //- Read a 2D table
        virtual void read
        (
            ISstream&,
            List<Tuple2<scalar, List<Tuple2<scalar, Type>>>>&
        ) const = 0;


public:

    //- Runtime type information
    TypeName("TableReader");

    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            TableReader,
            dictionary,
            (const dictionary& dict),
            (dict)
        );


    // Constructors

        //- Construct from dictionary
        TableReader(const dictionary& dict);

        //- Construct and return a clone
        virtual autoPtr<TableReader<Type>> clone() const = 0;


    // Selector
    static autoPtr<TableReader<Type>> New
    (
        const word& readerType,
        const dictionary& dict
    );


    //- Destructor
    virtual ~TableReader();


    // Member Functions

        //- Read a table
        template<class TableType>
        void operator()(const fileName&, TableType&) const;

        //- Write additional information
        virtual void write(Ostream& os) const;
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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "TableReader.C"
    #include "FoamTableReader.H"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
