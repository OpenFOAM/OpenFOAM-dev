/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2021-2026 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#ifndef compileTemplate_H
#define compileTemplate_H

#include "codedBase.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class compileTemplate Declaration
\*---------------------------------------------------------------------------*/

class compileTemplate
:
    public codedBase
{
    // Private Data

        //- Keywords associated with source code
        static const wordList codeKeys;

        //- Name of the dictionary variables in the source code
        static const wordList codeDictVars;

        //- Files to copy and filter and compile
        static wordList compileFiles_;

        //- List of template argument substitutions
        const List<Pair<word>> substitutions_;

        const dictionary dict_;


    // Private Member Functions

        fileName name(const word& instantiatedName) const;

        dictionary optionsDict
        (
            const word& templateName
        ) const;

        void setSubstitution(const Pair<word>& substitution);


public:

    // Constructors

        //- Construct from name and dictionary
        compileTemplate
        (
            const word& templateName,
            const word& instantiatedName,
            const List<Pair<word>>& substitutions
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
