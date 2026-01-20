/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2020-2026 OpenFOAM Foundation
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
    Foam::TableReaders::Embedded

Description
    Reads an interpolation table from the \c values entry in OpenFOAM-format.

    This is a list of Tuple2's where the first (x) column is scalar, and the
    second (y) column is the type to be interpolated.

Usage
    \verbatim
        values
        (
            (0.0 (1 2 3))
            (1.0 (4 5 6))
        );
    \endverbatim

SourceFiles
    Embedded.C

\*---------------------------------------------------------------------------*/

#ifndef EmbeddedTableReader_H
#define EmbeddedTableReader_H

#include "TableReader.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace TableReaders
{

/*---------------------------------------------------------------------------*\
                            Class Foam Declaration
\*---------------------------------------------------------------------------*/

template<class Coordinate, class Value>
class Embedded
:
    public TableReader<Coordinate, Value>
{
public:

    //- Runtime type information
    TypeName("embedded");


    // Constructors

        //- Default construct
        Embedded();

        //- Construct from name and dictionary
        Embedded
        (
            const word& name,
            const Function1s::unitConversions& units,
            const dictionary& dict
        );

        //- Construct and return a copy
        virtual autoPtr<TableReader<Coordinate, Value>> clone() const
        {
            return
                autoPtr<TableReader<Coordinate, Value>>
                (
                    new Embedded<Coordinate, Value>(*this)
                );
        }


    //- Destructor
    virtual ~Embedded();


    // Member Functions

        //- Read values
        virtual List<Tuple2<Coordinate, Value>> read
        (
            const Function1s::unitConversions& units,
            const dictionary& dict,
            const word& valuesKeyword="values"
        ) const;

        //- Read values
        List<Tuple2<Coordinate, Value>> read
        (
            const Function1s::unitConversions& units,
            Istream& is
        );

        //- Write settings and values
        virtual void write
        (
            Ostream& os,
            const Function1s::unitConversions& units,
            const List<Tuple2<Coordinate, Value>>& table,
            const word& valuesKeyword="values"
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace TableReaders
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "EmbeddedTableReader.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
