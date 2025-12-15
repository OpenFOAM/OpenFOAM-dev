/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2018-2025 OpenFOAM Foundation
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
    Foam::functionEntries::ifEntry

Description
    Conditional parsing of and within dictionary entries.

    E.g.
    \verbatim
        Uinlet 15;

        #if (#calc "${Uinlet} < 10")
            ..
        #else
            ..
        #endif

        pOutlet #if (#calc "${Uinlet} < 10") 0 #else 100 #endif;

        fluid       waterAndSludge;

        #if (#calc #{$fluid == "waterAndSludge" || $fluid == "waterAndSand" #})
            mixture     multicomponentMixture;
        #else
            mixture     pureMixture;
        #endif
    \endverbatim

    The condition should return any boolean representation supported by
    \c Switch i.e 0, 1, true, false, on, off, yes, no.

See also
    Foam::functionEntries::ifeqEntry
    Foam::functionEntries::elseEntry
    Foam::functionEntries::endifEntry
    Foam::Switch

SourceFiles
    ifEntry.C

\*---------------------------------------------------------------------------*/

#ifndef ifEntry_H
#define ifEntry_H

#include "ifeqEntry.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionEntries
{

/*---------------------------------------------------------------------------*\
                           Class ifEntry Declaration
\*---------------------------------------------------------------------------*/

class ifEntry
:
    public ifeqEntry
{
    friend class ifeqEntry;

    // Private Member Functions

        //- Execute the functionEntry in the given context
        template<class Context>
        bool execute
        (
            DynamicList<filePos>& stack,
            const dictionary& contextDict,
            Context& context,
            Istream& is
        ) const;


public:

    //- Runtime type information
    FunctionTypeName("#if");


    // Constructors

        //- Construct from line number, dictionary and Istream
        ifEntry
        (
            const label lineNumber,
            const dictionary& parentDict,
            Istream& is
        );

        //- Copy construct
        ifEntry(const ifEntry&) = default;

        //- Clone
        virtual autoPtr<entry> clone(const dictionary&) const
        {
            return autoPtr<entry>(new ifEntry(*this));
        }


    // Member Functions

        //- Expand the functionEntry into the contextDict
        virtual bool execute(dictionary& contextDict, Istream& is);

        //- Expand the functionEntry into the contextEntry
        static bool execute
        (
            const dictionary& contextDict,
            primitiveEntry& contextEntry,
            Istream&
        );


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const ifEntry&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionEntries
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "ifEntryTemplates.C"
    #include "ifeqEntryTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
