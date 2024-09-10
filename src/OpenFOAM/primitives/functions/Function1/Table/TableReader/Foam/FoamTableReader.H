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
    Foam::TableReaders::Foam

Description
    Reads an interpolation table from a file in OpenFOAM-format. This is a list
    of Tuples in which complex types are delimited by parentheses. By default,
    a Tuple2 is assumed in which the first column is the (scalar) x-axis and
    the second column is the y-axis. But finer control over a n-component Tuple
    can be achieved by specifying the column indices.

Usage
    A file containing the following will be read with the first element in each
    tuple as the x-coordinate and the second as the y-values:
    \verbatim
        (
            (0.0 (1 2 3))
            (1.0 (4 5 6))
        )
    \endverbatim

    A more complex table file might look as follows:
    \verbatim
        (
            (0 (1 2 3) 4 (5 6 7 8 9 10))
            (1 (2 3 4) 5 (6 7 8 9 10 11))
            (2 (3 4 5) 6 (7 8 9 10 11 12))
        )
    \endverbatim

    In this case, the columns to be used can be specified using the \c columns
    entry. In this example, the second scalar column is taken as the x-axis,
    and the six-component symmTensor-s are used as the values:
    \verbatim
        columns     (2 3);
    \endverbatim

SourceFiles
    FoamTableReader.C

\*---------------------------------------------------------------------------*/

#ifndef FoamTableReader_H
#define FoamTableReader_H

#include "TableFileReader.H"
#include "labelPair.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace TableReaders
{

/*---------------------------------------------------------------------------*\
                            Class Foam Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class Foam
:
    public TableFileReader<Type>
{
    // Private Data

        //- Column indices
        const labelPair columns_;


    // Private Member Functions

        //- Read a 1D table
        virtual void read(ISstream&, List<Tuple2<scalar, Type>> &) const;


public:

    //- Runtime type information
    TypeName("foam");


    // Constructors

        //- Construct from name and dictionary
        Foam
        (
            const word& name,
            const Function1s::unitConversions& units,
            const dictionary& dict
        );

        //- Construct and return a copy
        virtual autoPtr<TableReader<Type>> clone() const
        {
            return autoPtr<TableReader<Type>>(new Foam<Type>(*this));
        }


    //- Destructor
    virtual ~Foam();


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
    #include "FoamTableReader.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
