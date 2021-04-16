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

#include "basicThermo.H"
#include "compileThermo.H"
#include "dynamicCodeContext.H"
#include "Time.H"
#include "IFstream.H"
#include "OSspecific.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

template<>
const Foam::wordList Foam::CodedBase<Foam::compileThermo>::codeKeys_ =
{};


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::dictionary Foam::compileThermo::optionsDict
(
    const word& thermoTypeName,
    const word& instantiatedThermoName
) const
{
    IFstream optionsFile(dynamicCode::resolveTemplate(thermoTypeName));
    if (!optionsFile.good())
    {
        FatalErrorInFunction
            << "Failed to open dictionary file " << thermoTypeName
            << exit(FatalError);
    }

    dictionary dict(optionsFile);

    fileName thermoTypeFileName(instantiatedThermoName);
    thermoTypeFileName.replaceAll(',', '_');
    thermoTypeFileName.replaceAll('<', '_');
    thermoTypeFileName.replaceAll('>', '_');

    dict.add("name", thermoTypeFileName);

    return dict;
}


void Foam::compileThermo::setFilterVariable
(
    dynamicCode& dynCode,
    const dynamicCodeContext& context,
    const word& name,
    const word& type
) const
{
    const HashSet<word> types(context.dict().lookup(name));
    if (!types.found(type))
    {
        FatalIOErrorInFunction(thermoTypeDict_)
            << "Unknown " << name << " type " << type << nl
            << "Supported " << name << " types: " << types
            << exit(FatalIOError);
    }

    dynCode.setFilterVariable(name, type);
}


void Foam::compileThermo::setFilterVariable
(
    dynamicCode& dynCode,
    const dynamicCodeContext& context,
    const word& name
) const
{
    setFilterVariable
    (
        dynCode,
        context,
        name,
        thermoTypeDict_.lookup<word>(name)
    );
}


void Foam::compileThermo::setFilterRenamedVariable
(
    dynamicCode& dynCode,
    const dynamicCodeContext& context,
    const word& name,
    const word& typeRenameMapName
) const
{
    word type(thermoTypeDict_.lookup<word>(name));

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

    setFilterVariable(dynCode, context, name, type);
}


void Foam::compileThermo::prepare
(
    dynamicCode& dynCode,
    const dynamicCodeContext& context
) const
{
    dynCode.setFilterVariable("typeName", codeName());

    const word type(thermoTypeDict_.lookup<word>("type"));
    dynCode.setFilterVariable("type", type);

    const HashTable<word> typeToBaseMap(context.dict().lookup("baseType"));
    dynCode.setFilterVariable("baseType", typeToBaseMap[type]);

    setFilterRenamedVariable(dynCode, context, "mixture", "renameMixture");
    setFilterVariable(dynCode, context, "transport");
    setFilterVariable(dynCode, context, "thermo");
    setFilterVariable(dynCode, context, "equationOfState");
    setFilterVariable(dynCode, context, "energy");

    // Compile filtered C template
    dynCode.addCompileFile(thermoTypeName_ + ".C");

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

Foam::compileThermo::compileThermo
(
    const word& thermoTypeName,
    const word& instantiatedThermoName,
    const dictionary& thermoTypeDict
)
:
    CodedBase<compileThermo>
    (
        optionsDict(thermoTypeName, instantiatedThermoName)
    ),
    thermoTypeName_(thermoTypeName),
    thermoTypeDict_(thermoTypeDict)
{
    this->updateLibrary();
}


// ************************************************************************* //
