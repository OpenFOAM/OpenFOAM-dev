/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2024 OpenFOAM Foundation
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
    Foam::namedUnitConversion

Description
    ...

SourceFiles
    ...

\*---------------------------------------------------------------------------*/

#ifndef namedUnitConversion_H
#define namedUnitConversion_H

#include "unitConversion.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class Istream;
class Ostream;
class entry;

// Forward declaration of friend functions and operators
class namedUnitConversion;
Istream& operator>>(Istream&, namedUnitConversion&);
Ostream& operator<<(Ostream&, const namedUnitConversion&);

/*---------------------------------------------------------------------------*\
                          Class namedUnitConversion Declaration
\*---------------------------------------------------------------------------*/

class namedUnitConversion
:
    public unitConversion
{
private:

    // Private Data

        //- The name
        string name_;


    // Private Member Functions

        //- Read the name
        void readName(const entry& e);

        //- Read the name
        void readName(const word& keyword, const dictionary&);


public:

    // Static Data Members

        //- Run time type information
        ClassName("namedUnitConversion");


    // Constructors

        //- Construct from name
        namedUnitConversion(const string&);

        //- Copy constructor
        namedUnitConversion(const namedUnitConversion&) = default;

        //- Move constructor
        namedUnitConversion(namedUnitConversion&&) = default;


    // Member Functions

        //- Update
        void read(const word& keyword, const dictionary&);

        //- Update if found in the dictionary
        bool readIfPresent(const word& keyword, const dictionary&);


    // IOstream Operators

        //- Read from stream
        friend Istream& operator>>(Istream&, namedUnitConversion&);

        //- Write to stream
        friend Ostream& operator<<(Ostream&, const namedUnitConversion&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
