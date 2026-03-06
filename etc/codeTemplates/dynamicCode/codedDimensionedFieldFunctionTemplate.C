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

\*---------------------------------------------------------------------------*/

#include "codedDimensionedFieldFunctionTemplate.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace DimensionedFieldFunctions
{
    defineTypeNameAndDebug
    (
        ${typeName}DimensionedFieldFunction${DimensionedFieldTypeName},
        0
    );
}

DimensionedFieldFunction<${DimensionedFieldType}>::
adddictionaryConstructorToTable<DimensionedFieldFunctions::
${typeName}DimensionedFieldFunction${DimensionedFieldTypeName}>
${typeName}DimensionedFieldFunction${DimensionedFieldTypeName}ConstructorToTable_;

}


// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

extern "C"
{
    // Unique function name that can be checked
    // to ensure the correct library version has been loaded
    void ${uniqueFunctionName}(bool load)
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


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::DimensionedFieldFunctions::
${typeName}DimensionedFieldFunction${DimensionedFieldTypeName}::
${typeName}DimensionedFieldFunction${DimensionedFieldTypeName}
(
    const dictionary& dict,
    ${DimensionedFieldType}& field_
)
:
    DimensionedFieldFunction<${DimensionedFieldType}>(dict, field_),
    field(field_)
{
    if (${verbose})
    {
        Info<< "Construct ${typeName} sha1: ${SHA1sum} from dictionary\n";
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::DimensionedFieldFunctions::
${typeName}DimensionedFieldFunction${DimensionedFieldTypeName}::
~${typeName}DimensionedFieldFunction${DimensionedFieldTypeName}()
{
    if (${verbose})
    {
        Info<< "Destroy ${typeName} sha1: ${SHA1sum}\n";
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::DimensionedFieldFunctions::
${typeName}DimensionedFieldFunction${DimensionedFieldTypeName}::evaluate()
{
    const DimensionedField<vector, GeoMesh, Field>& C(field.mesh().C());

//{{{ begin code
    ${evaluate}
//}}} end code
}


void Foam::DimensionedFieldFunctions::
${typeName}DimensionedFieldFunction${DimensionedFieldTypeName}::update()
{
//{{{ begin code
    ${update}
//}}} end code
}


void Foam::DimensionedFieldFunctions::
${typeName}DimensionedFieldFunction${DimensionedFieldTypeName}::
write
(
    Ostream& os
) const
{
    NotImplemented;
}


// ************************************************************************* i/
