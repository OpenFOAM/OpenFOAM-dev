/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2020 OpenFOAM Foundation
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

    Usage:
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

template<class Type>
class Embedded
:
    public TableReader<Type>
{

public:

    //- Runtime type information
    TypeName("embedded");


    // Constructors

        //- Construct from dictionary
        Embedded
        (
            const word& name,
            const dictionary &dict,
            List<Tuple2<scalar, Type>>& table
        );

        //- Construct and return a copy
        virtual autoPtr<TableReader<Type>> clone() const
        {
            return autoPtr<TableReader<Type>>(new Embedded<Type>(*this));
        }


    //- Destructor
    virtual ~Embedded();


    // Member Functions

        //- Write the table values
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
    #include "EmbeddedTableReader.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
