/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012-2013 OpenFOAM Foundation
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
    Foam::doxygenXmlParser

Description
    Parser for doxygen XML

SourceFiles
    doxygenXmlParser.C

\*---------------------------------------------------------------------------*/

#ifndef doxygenXmlParser_H
#define doxygenXmlParser_H

#include "dictionary.H"
#include "IFstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class doxygenXmlParser Declaration
\*---------------------------------------------------------------------------*/

class doxygenXmlParser
:
    public dictionary
{

public:

    //- Construct from components
    doxygenXmlParser
    (
        const fileName& fName,
        const string& startTag,
        const string& searchStr,
        const bool exactMatch
    );


    // Member functions

        //- Skip past a block
        void skipBlock(IFstream& is, const word& blockName) const;

        //- Skip forward to block
        void skipForward(IFstream& is, const word& blockName) const;

        //- Return the entry
        template<class Type>
        void getEntry(IFstream& is, Type& entry) const;

        //- Return the entry value
        template<class Type>
        void getValue(IFstream& is, Type& entry) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "doxygenXmlParserTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
