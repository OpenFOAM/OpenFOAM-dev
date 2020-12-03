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
    Template for use with dynamic code generation of a Function2.

    - without state

SourceFiles
    Function2Template.C

\*---------------------------------------------------------------------------*/

#ifndef Function2Template_H
#define Function2Template_H

#include "Function2.H"

//{{{ begin codeInclude
${codeInclude}
//}}} end codeInclude

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace Function2s
{

/*---------------------------------------------------------------------------*\
                           A templated Function2
\*---------------------------------------------------------------------------*/

class ${typeName}Function2${TemplateType}
:
    public FieldFunction2<${TemplateType}, ${typeName}Function2${TemplateType}>
{

public:

    // Runtime type information
    TypeName("${typeName}");


    // Constructors

        //- Construct from entry name and dictionary
        ${typeName}Function2${TemplateType}
        (
            const word& entryName,
            const dictionary& dict
        );

        //- Copy constructor
        ${typeName}Function2${TemplateType}
        (
            const ${typeName}Function2${TemplateType}& f1
        );

        //- Construct and return a clone
        virtual tmp<Function2<${TemplateType}>> clone() const
        {
            return tmp<Function2<${TemplateType}>>
            (
                new ${typeName}Function2${TemplateType}(*this)
            );
        }


    //- Destructor
    virtual ~${typeName}Function2${TemplateType}();


    // Member Functions

        //- Return value as a function of two scalars
        inline virtual ${TemplateType} value
        (
            const scalar x,
            const scalar y
        ) const
        {
//{{{ begin code
    ${code}
//}}} end code
        }

        //- Write data to dictionary stream
        virtual void write(Ostream& os) const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const ${typeName}Function2${TemplateType}&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Function2s
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
