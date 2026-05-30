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

Description
    Template for use with dynamic code generation of a zoneGenerator.

SourceFiles
    zoneGeneratorTemplate.C

\*---------------------------------------------------------------------------*/

#ifndef codedZoneGeneratorTemplate_H
#define codedZoneGeneratorTemplate_H

#include "codedZoneGenerator.H"

//{{{ begin codeInclude
${codeInclude}
//}}} end codeInclude

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                              A Coded zoneGenerator
\*---------------------------------------------------------------------------*/

class ${typeName}ZoneGenerator
:
    public zoneGenerator
{
    // Private Data


public:

    //- Runtime type information
    TypeName("${typeName}");


    // Constructors

        //- Construct from Time and dictionary
        ${typeName}ZoneGenerator
        (
            const word& name,
            const polyMesh& mesh,
            const dictionary&
        );

        //- Disallow default bitwise copy construction
        ${typeName}ZoneGenerator(const ${typeName}ZoneGenerator&) = delete;


    //- Destructor
    virtual ~${typeName}ZoneGenerator();


    // Member Functions

        //- Generate and return the zoneSet
        virtual zoneSet generate() const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const ${typeName}ZoneGenerator&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
