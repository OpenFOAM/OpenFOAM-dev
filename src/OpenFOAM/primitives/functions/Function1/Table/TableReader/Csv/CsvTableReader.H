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
    Foam::TableReaders::Csv

Description
    Reads an interpolation table from a file in CSV-format. Entries govern the
    layout of the CSV file. The index of the first (x) column of the table is
    given by the refColumn entry, and is always scalar. The indices of the
    components of the second (y) column are given by the componentColumns
    entry.

    Usage:
    \verbatim
        nHeaderLine         4;          // number of header lines
        refColumn           0;          // reference column index
        componentColumns    (1 2 3);    // component column indices
        separator           ",";        // optional (defaults to ",")
        mergeSeparators     no;         // merge multiple separators
    \endverbatim

SourceFiles
    CsvTableReader.C

\*---------------------------------------------------------------------------*/

#ifndef CsvTableReader_H
#define CsvTableReader_H

#include "TableFileReader.H"
#include "labelList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace TableReaders
{

/*---------------------------------------------------------------------------*\
                             Class Csv Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class Csv
:
    public TableFileReader<Type>
{
    // Private Data

        //- Number of header lines
        const label nHeaderLine_;

        //- Reference column
        const label refColumn_;

        //- Labels of the components
        const labelList componentColumns_;

        //- Separator character
        const char separator_;

        //- Merge separators flag; e.g. ',,,' becomes ','
        bool mergeSeparators_;


    // Private Member functions

        //- Read the next value from the split string
        Type readValue(const List<string>&) const;

        //- Read a 1D table
        virtual void read(ISstream&, List<Tuple2<scalar, Type>>&) const;


public:

    //- Runtime type information
    TypeName("csv");


    // Constructors

        //- Construct from dictionary
        Csv
        (
            const word& name,
            const dictionary &dict,
            List<Tuple2<scalar, Type>>& table
        );

        //- Construct and return a copy
        virtual autoPtr<TableReader<Type>> clone() const
        {
            return autoPtr<TableReader<Type>>(new Csv<Type>(*this));
        }


    //- Destructor
    virtual ~Csv();


    // Member Functions

        //- Write the CSV parameters
        virtual void write
        (
            Ostream& os,
            const List<Tuple2<scalar, Type>>& table
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace TableReaders
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "CsvTableReader.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
