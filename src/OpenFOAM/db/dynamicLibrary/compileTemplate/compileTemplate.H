/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2021 OpenFOAM Foundation
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
#include "dynamicCode.H"

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
    // Private Member Data

        //- Name of the class Thermo basicThermo is instantiated on
        const word templateName_;

        //- List of template argument substitutions
        const List<Pair<word>> substitutions_;


    // Private Member Functions

        fileName name(const word& instantiatedName) const;

        dictionary optionsDict
        (
            const word& templateName
        ) const;

        void setFilterVariable
        (
            dynamicCode& dynCode,
            const dynamicCodeContext& context,
            const Pair<word>& substitution
        ) const;

        //- Get the keywords associated with source code
        virtual wordList codeKeys() const
        {
            return wordList::null();
        }

        //- Adapt the context for the current object
        virtual void prepare(dynamicCode&, const dynamicCodeContext&) const;

        //- Clear the ptr to the redirected object
        virtual void clearRedirect() const
        {}


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
