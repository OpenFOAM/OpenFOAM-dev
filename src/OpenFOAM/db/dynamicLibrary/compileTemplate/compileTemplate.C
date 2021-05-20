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

#include "compileTemplate.H"
#include "dynamicCodeContext.H"
#include "Time.H"
#include "IFstream.H"
#include "OSspecific.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::fileName Foam::compileTemplate::name
(
    const word& instantiatedName
) const
{
    fileName templateFileName(instantiatedName);
    templateFileName.replaceAll(',', '_');
    templateFileName.replaceAll('<', '_');
    templateFileName.replaceAll('>', '_');

    return templateFileName;
}


Foam::dictionary Foam::compileTemplate::optionsDict
(
    const word& templateName
) const
{
    IFstream optionsFile(dynamicCode::resolveTemplate(templateName));
    if (!optionsFile.good())
    {
        FatalErrorInFunction
            << "Failed to open dictionary file " << templateName
            << exit(FatalError);
    }

    return dictionary(optionsFile);
}


void Foam::compileTemplate::setFilterVariable
(
    dynamicCode& dynCode,
    const dynamicCodeContext& context,
    const Pair<word>& substitution
) const
{
    const word& name(substitution.first());
    word type(substitution.second());
    const word typeRenameMapName(name + "Renamed");

    if (context.dict().found(name))
    {
        const HashSet<word> types(context.dict().lookup(name));
        if (!types.found(type))
        {
            FatalIOErrorInFunction(context.dict())
                << "Unknown " << name << " type " << type << nl
                << "Supported " << name << " types: " << types
                << exit(FatalIOError);
        }
    }

    if (context.dict().found(typeRenameMapName))
    {
        const HashTable<word> renameMap
        (
            context.dict().lookup(typeRenameMapName)
        );

        if (renameMap.found(type))
        {
            type = renameMap[type];
        }
    }

    dynCode.setFilterVariable(name, type);

    const word typeBase(name + "Base");
    if (context.dict().found(typeBase))
    {
        const HashTable<word> typeToBaseMap(context.dict().lookup(typeBase));
        dynCode.setFilterVariable(typeBase, typeToBaseMap[type]);
    }
}


void Foam::compileTemplate::prepare
(
    dynamicCode& dynCode,
    const dynamicCodeContext& context
) const
{
    dynCode.setFilterVariable("typeName", codeName());

    forAll(substitutions_, i)
    {
        setFilterVariable(dynCode, context, substitutions_[i]);
    }

    // Compile filtered C template
    dynCode.addCompileFile(codeTemplateC(templateName_));

    // Define Make/options
    dynCode.setMakeOptions(context.options() + "\n\n" + context.libs());

    // Debugging: make verbose
    if (debug)
    {
        dynCode.setFilterVariable("verbose", "true");
        Info<<"compile " << codeName() << " sha1: "
            << context.sha1() << endl;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::compileTemplate::compileTemplate
(
    const word& templateName,
    const word& instantiatedName,
    const List<Pair<word>>& substitutions
)
:
    codedBase(name(instantiatedName), optionsDict(templateName)),
    templateName_(templateName),
    substitutions_(substitutions)
{
    this->updateLibrary();
}


// ************************************************************************* //
