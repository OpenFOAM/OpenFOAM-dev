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

#include "Function1Template.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace Function1s
{
    defineTypeNameAndDebug(${typeName}Function1${TemplateType}, 0);
}
    Function1<${TemplateType}>::adddictionaryConstructorToTable<Function1s::
        ${typeName}Function1${TemplateType}>
        ${typeName}Function1${TemplateType}ConstructorToTable_;
}


// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

extern "C"
{
    // dynamicCode:
    // SHA1 = ${SHA1sum}
    //
    // Unique function name that can be checked if the correct library version
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


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::Function1s::${typeName}Function1${TemplateType}::
${typeName}Function1${TemplateType}
(
    const word& entryName,
    const dictionary& dict
)
:
    FieldFunction1<${TemplateType}, ${typeName}Function1${TemplateType}>
    (
        entryName
    )
{
    if (${verbose:-false})
    {
        Info<< "Construct ${typeName} sha1: ${SHA1sum} from dictionary\n";
    }
}


Foam::Function1s::${typeName}Function1${TemplateType}::
${typeName}Function1${TemplateType}
(
    const ${typeName}Function1${TemplateType}& f1
)
:
    FieldFunction1<${TemplateType}, ${typeName}Function1${TemplateType}>
    (
        f1
    )
{
    if (${verbose:-false})
    {
        Info<< "Construct ${typeName} sha1: ${SHA1sum} as copy\n";
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::Function1s::${typeName}Function1${TemplateType}::
~${typeName}Function1${TemplateType}()
{
    if (${verbose:-false})
    {
        Info<< "Destroy ${typeName} sha1: ${SHA1sum}\n";
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::${TemplateType}
Foam::Function1s::${typeName}Function1${TemplateType}::integral
(
    const scalar x1,
    const scalar x2
) const
{
    NotImplemented;
    return pTraits<${TemplateType}>::zero;
}


void Foam::Function1s::${typeName}Function1${TemplateType}::write
(
    Ostream& os
) const
{
    NotImplemented;
}


// ************************************************************************* i/
