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
    Foam::TableReaders::Foam

Description
    Reads an interpolation table from a file in OpenFOAM-format. This is a list
    of Tuple2's where the first (x) column is scalar, and the second (y) column
    is the type to be interpolated.

    Usage:
    \verbatim
        (
            (0.0 (1 2 3))
            (1.0 (4 5 6))
        );
    \endverbatim

SourceFiles
    FoamTableReader.C

\*---------------------------------------------------------------------------*/

#ifndef FoamTableReader_H
#define FoamTableReader_H

#include "TableFileReader.H"

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
    // Private Member Functions

        //- Read a 1D table
        virtual void read(ISstream&, List<Tuple2<scalar, Type>> &) const;


public:

    //- Runtime type information
    TypeName("foam");


    // Constructors

        //- Construct from dictionary
        Foam
        (
            const word& name,
            const dictionary &dict,
            List<Tuple2<scalar, Type>>& table
        );

        //- Construct and return a copy
        virtual autoPtr<TableReader<Type>> clone() const
        {
            return autoPtr<TableReader<Type>>(new Foam<Type>(*this));
        }


    //- Destructor
    virtual ~Foam();
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
