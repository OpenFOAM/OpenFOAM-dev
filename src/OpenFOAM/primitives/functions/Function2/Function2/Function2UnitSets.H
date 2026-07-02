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
    Foam::Function2

Description
    Struct containing three unitSets for use in converting both the two
    arguments and the value of a Function2.

SourceFiles
    Function2UnitSets.C

\*---------------------------------------------------------------------------*/

#ifndef Function2UnitSets_H
#define Function2UnitSets_H

#include "unitSet.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace Function2s
{

/*---------------------------------------------------------------------------*\
                 Class Function2s::unitSets Declaration
\*---------------------------------------------------------------------------*/

struct unitSets
{
    // Public Data

        //- Unit conversion for x-axis values
        unitSet x;

        //- Unit conversion for y-axis values
        unitSet y;

        //- Unit conversion for result values
        unitSet value;


    // Constructors

        //- Construct from initialiser list
        unitSets(std::initializer_list<unitSet>);

        //- Construct from stream
        unitSets(Istream& is);


    // Member Functions

        //- Update from stream
        void read(Istream& is);
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Function2s

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Global Functions

//- Generate an error in an context where unit conversions are not supported
void assertNoConvertUnits
(
    const word& typeName,
    const Function2s::unitSets& units,
    const dictionary& dict
);

//- Write a units entry
void writeEntry(Ostream& os, const Function2s::unitSets& units);

//- Write a type with Function2 unit conversion
template<class Type>
inline void writeEntry(Ostream& os, const Function2s::unitSets&, const Type&);


// IOStream Operators

//- Read from stream
Istream& operator>>(Istream& is, Function2s::unitSets&);

//- Write to stream
Ostream& operator<<(Ostream& os, const Function2s::unitSets&);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "Function2UnitSetsTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
