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
    Foam::TableReaders::Csv

Description
    Reads an interpolation table from a file in CSV-format. Entries govern the
    layout of the CSV file. The indices of the columns of the table that are to
    be used are given by the columns entry. This is a tuple for which the first
    part is the index of the column used for the x-axis, and the second part is
    the column index used for the scalar values, or a list of column indices
    used for the components of vector, tensor, etc..., values.

Usage
    \verbatim
        nHeaderLine         4;           // number of header lines
        columns             (0 (1 2 3)); // column indices for vector values
        separator           ",";         // optional (defaults to ",")
        mergeSeparators     no;          // merge multiple separators
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
                           Using CsvVoid Declaration
\*---------------------------------------------------------------------------*/

template<typename T>
using CsvVoid = void;


/*---------------------------------------------------------------------------*\
                        Struct CsvLabelType Declaration
\*---------------------------------------------------------------------------*/

template<class Type, typename = void>
struct CsvLabelType
{
    typedef label type;

    typedef FixedList<label, 1> oldType;

    type operator()(const oldType& ot)
    {
        return ot[0];
    }
};


template<class Type>
struct CsvLabelType<Type, CsvVoid<typename pTraits<Type>::labelType>>
{
    typedef typename pTraits<Type>::labelType type;

    typedef FixedList<label, pTraits<Type>::nComponents> oldType;

    type operator()(const oldType& ot) const
    {
        type t;
        for (direction i = 0; i < pTraits<Type>::nComponents; ++ i)
        {
            setComponent(t, i) = ot[i];
        }
        return t;
    }
};


/*---------------------------------------------------------------------------*\
                             Class Csv Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class Csv
:
    public TableFileReader<Type>
{
    // Private Typedefs

        //- Type of the column indices
        typedef Tuple2<label, typename CsvLabelType<Type>::type> columnIndices;


    // Private Data

        //- Number of header lines
        const label nHeaderLine_;

        //- Column indices
        const columnIndices columns_;

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

        //- Construct from name and dictionary
        Csv
        (
            const word& name,
            const Function1s::unitConversions& units,
            const dictionary& dict
        );

        //- Construct and return a copy
        virtual autoPtr<TableReader<Type>> clone() const
        {
            return autoPtr<TableReader<Type>>(new Csv<Type>(*this));
        }


    //- Destructor
    virtual ~Csv();


    // Member Functions

        //- Write settings and values
        virtual void write
        (
            Ostream& os,
            const Function1s::unitConversions& units,
            const List<Tuple2<scalar, Type>>& table,
            const word& valuesKeyword=word::null
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
