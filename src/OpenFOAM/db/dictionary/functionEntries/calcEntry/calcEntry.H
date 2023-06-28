/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2023 OpenFOAM Foundation
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
    Foam::functionEntries::calcEntry

Description
    Compiles and executes code string expressions,
    returning the result to the dictionary

Usage
    \verbatim
    a 1.1;
    b 3.2;
    c #calc "$a*$b";
    \endverbatim

    Special care is required for calc entries that include a division since
    "/" is also used as the scoping operator to identify keywords in
    sub-dictionaries. For example, "$a/b" expects a keyword "b" within a
    sub-dictionary named "a". A division can be correctly executed by using a
    space between a variables and "/", e.g.

    \verbatim
    c #calc "$a / $b";
    \endverbatim

    or bracketing the variable, e.g.

    \verbatim
    c #calc "($a)/$b";
    \endverbatim

    The code string can be delimited either by '"..."' in which newlines must be
    escaped or '#{...#}' which directly supports multi-line strings and is more
    convenient when evaluating string expressions by avoiding the need to
    escape the quotes, e.g.

    \verbatim
    c #calc #{ $a*sqrt($b) #};
    \endverbatim

    \verbatim
    s "field";
    fieldName #calc #{ $<string>s + "Name" #};
    \endverbatim

    Additional include files for the #calc code compilation can be specified
    using the #calcInclude entry, e.g. if functions from transform.H are used
    \verbatim
    angleOfAttack   5; // degs

    angle           #calc "-degToRad($angleOfAttack)";

    #calcInclude    "transform.H"
    liftDir         #calc "transform(Ry($angle), vector(0, 0, 1))";
    dragDir         #calc "transform(Ry($angle), vector(1, 0, 0))";
    \endverbatim

    Note:
        Internally this is just a wrapper around codeStream functionality - the
        #calc string is used to construct a dictionary for codeStream.

See also
    Foam::functionEntries::calcIncludeEntry

SourceFiles
    calcEntry.C

\*---------------------------------------------------------------------------*/

#ifndef calcEntry_H
#define calcEntry_H

#include "functionEntry.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class dlLibraryTable;

namespace functionEntries
{

/*---------------------------------------------------------------------------*\
                          Class calcEntry Declaration
\*---------------------------------------------------------------------------*/

class calcEntry
:
    public functionEntry
{
    // Private Member Functions

        //- Perform the calculation and return the resulting string
        static string calc
        (
            const dictionary& dict,
            Istream& is
        );


public:

    //- Runtime type information
    ClassName("calc");


    // Constructors

        //- Disallow default bitwise copy construction
        calcEntry(const calcEntry&) = delete;


    // Member Functions

        //- Execute the functionEntry in a sub-dict context
        static bool execute(dictionary& dict, Istream&);

        //- Execute the functionEntry in a primitiveEntry context
        static bool execute
        (
            const dictionary& dict,
            primitiveEntry&,
            Istream&
        );


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const calcEntry&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionEntries
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
