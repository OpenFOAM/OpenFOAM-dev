/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2023 OpenFOAM Foundation
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
    Foam::functionEntries::calcIncludeEntry

Description
    Specify an include file for #calc, expects a single string to follow.

    For example if functions from transform.H are used in the #calc expression
    \verbatim
    angleOfAttack   5; // degs

    angle           #calc "-degToRad($angleOfAttack)";

    #calcInclude    "transform.H"
    liftDir         #calc "transform(Ry($angle), vector(0, 0, 1))";
    dragDir         #calc "transform(Ry($angle), vector(1, 0, 0))";
    \endverbatim

    The usual expansion of environment variables and other constructs
    (eg, the \c ~OpenFOAM/ expansion) is retained.

See also
    Foam::functionEntries::calcEntry

SourceFiles
    calcIncludeEntry.C

\*---------------------------------------------------------------------------*/

#ifndef calcIncludeEntry_H
#define calcIncludeEntry_H

#include "functionEntry.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionEntries
{

/*---------------------------------------------------------------------------*\
                        Class calcIncludeEntry Declaration
\*---------------------------------------------------------------------------*/

class calcIncludeEntry
:
    public functionEntry
{
    // Private static data

        // Static include file name cache
        static DynamicList<fileName> includeFiles_;


public:

    //- Runtime type information
    ClassName("calcInclude");


    // Constructors

        //- Disallow default bitwise copy construction
        calcIncludeEntry(const calcIncludeEntry&) = delete;


    // Member Functions

        //- Execute the functionEntry in a sub-dict context
        static bool execute(dictionary& parentDict, Istream&);

        //- Reset the cache of #calc include file names
        static void clear();

        //- Add the cached include file names to the codeInclude entry
        //  in codeDict
        static void codeInclude(dictionary& codeDict);


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const calcIncludeEntry&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionEntries
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
