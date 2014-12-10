/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013-2014 OpenFOAM Foundation
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

#include "codedFvOptionTemplate.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "unitConversion.H"
#include "fvMatrix.H"

//{{{ begin codeInclude
${codeInclude}
//}}} end codeInclude


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace fv
{

// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

//{{{ begin localCode
${localCode}
//}}} end localCode


// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

extern "C"
{
    // dynamicCode:
    // SHA1 = ${SHA1sum}
    //
    // unique function name that can be checked if the correct library version
    // has been loaded
    void ${typeName}_${SHA1sum}(bool load)
    {
        if (load)
        {
            // code that can be explicitly executed after loading
        }
        else
        {
            // code that can be explicitly executed before unloading
        }
    }
}

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

//makeRemovablePatchTypeField
//(
//    fvPatch${FieldType},
//    ${typeName}FvOption${SourceType}
//);
defineTypeNameAndDebug(${typeName}FvOption${SourceType}, 0);
addRemovableToRunTimeSelectionTable
(
    option,
    ${typeName}FvOption${SourceType},
    dictionary
);


const char* const ${typeName}FvOption${SourceType}::SHA1sum =
    "${SHA1sum}";


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

${typeName}FvOption${SourceType}::
${typeName}FvOption${SourceType}
(
    const word& name,
    const word& modelType,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    option(name, modelType, dict, mesh)
{
    if (${verbose:-false})
    {
        Info<<"construct ${typeName} sha1: ${SHA1sum}"
            " from components\n";
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

${typeName}FvOption${SourceType}::
~${typeName}FvOption${SourceType}()
{
    if (${verbose:-false})
    {
        Info<<"destroy ${typeName} sha1: ${SHA1sum}\n";
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void ${typeName}FvOption${SourceType}::correct
(
    GeometricField<${TemplateType}, fvPatchField, volMesh>& fld
)
{
    if (${verbose:-false})
    {
        Info<<"${typeName}FvOption${SourceType}::correct()\n";
    }

//{{{ begin code
    ${codeCorrect}
//}}} end code
}


void ${typeName}FvOption${SourceType}::addSup
(
    fvMatrix<${TemplateType}>& eqn,
    const label fieldI
)
{
    if (${verbose:-false})
    {
        Info<<"${typeName}FvOption${SourceType}::addSup()\n";
    }

//{{{ begin code
    ${codeAddSup}
//}}} end code
}


void ${typeName}FvOption${SourceType}::addSup
(
    const volScalarField& rho,
    fvMatrix<${TemplateType}>& eqn,
    const label fieldI
)
{
    if (${verbose:-false})
    {
        Info<<"${typeName}FvOption${SourceType}::addSup()\n";
    }

//{{{ begin code
    ${codeAddSup}
//}}} end code
}


void ${typeName}FvOption${SourceType}::setValue
(
    fvMatrix<${TemplateType}>& eqn,
    const label fieldI
)
{
    if (${verbose:-false})
    {
        Info<<"${typeName}FvOption${SourceType}::setValue()\n";
    }

//{{{ begin code
    ${codeSetValue}
//}}} end code
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

} // End namespace fv
// ************************************************************************* //
