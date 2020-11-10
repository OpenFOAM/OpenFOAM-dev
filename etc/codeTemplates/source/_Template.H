/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) YEAR OpenFOAM Foundation
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
    Foam::CLASSNAME

Description

SourceFiles
    CLASSNAMEI.H
    CLASSNAME.C
    CLASSNAMEIO.C

\*---------------------------------------------------------------------------*/

#ifndef CLASSNAME_H
#define CLASSNAME_H

#include ".H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class Istream;
class Ostream;

// Forward declaration of friend functions and operators
class CLASSNAME;
Istream& operator>>(Istream&, CLASSNAME&);
Ostream& operator<<(Ostream&, const CLASSNAME&);


/*---------------------------------------------------------------------------*\
                          Class CLASSNAME Declaration
\*---------------------------------------------------------------------------*/

class CLASSNAME
:
    public baseClassName
{
    // Private Data

        //- Description of data_
        dataType data_;


    // Private Member Functions


public:

    // Static Data Members

        //- Static data staticData
        static const dataType staticData;


    // Constructors

        //- Construct null
        CLASSNAME();

        //- Construct from components
        CLASSNAME(const dataType& data);

        //- Construct from Istream
        CLASSNAME(Istream&);

        //- Copy constructor
        CLASSNAME(const CLASSNAME&);

        //- Disallow default bitwise copy construction
        CLASSNAME(const CLASSNAME&) = delete;


    // Selectors

        //- Select null constructed
        static autoPtr<CLASSNAME> New();


    //- Destructor
    ~CLASSNAME();


    // Member Functions

        // Access

        // Check

        // Edit

        // Write


    // Member Operators

        void operator=(const CLASSNAME&);

        //- Disallow default bitwise assignment
        void operator=(const CLASSNAME&) = delete;


    // Friend Functions

    // Friend Operators

    // IOstream Operators

        friend Istream& operator>>(Istream&, CLASSNAME&);
        friend Ostream& operator<<(Ostream&, const CLASSNAME&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "CLASSNAMEI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
