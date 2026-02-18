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

#include "compileTemplate.H"
#include "Time.H"
#include "IFstream.H"
#include "OSspecific.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const Foam::wordList Foam::compileTemplate::codeKeys(wordList::null());

const Foam::wordList Foam::compileTemplate::codeDictVars(wordList::null());

Foam::wordList Foam::compileTemplate::compileFiles_(wordList::null());


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
    const Pair<word>& substitution
) const
{
    const word& name(substitution.first());
    word type(substitution.second());
    const word typeRenameMapName(name + "Renamed");

    if (dict_.found(name))
    {
        const HashSet<word> types(dict_.lookup(name));
        if (!types.found(type))
        {
            FatalIOErrorInFunction(dict_)
                << "Unknown " << name << " type " << type << nl
                << "Supported " << name << " types: " << types
                << exit(FatalIOError);
        }
    }

    if (dict_.found(typeRenameMapName))
    {
        const HashTable<word> renameMap
        (
            dict_.lookup(typeRenameMapName)
        );

        if (renameMap.found(type))
        {
            type = renameMap[type];
        }
    }

    codedBase::setFilterVariable(name, type);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::compileTemplate::compileTemplate
(
    const word& templateName,
    const word& instantiatedName,
    const List<Pair<word>>& substitutions
)
:
    codedBase
    (
        name(instantiatedName),
        optionsDict(templateName),
        codeKeys,
        codeDictVars,
        word::null,
        compileFiles_,
        wordList::null()
    ),
    substitutions_(substitutions),
    dict_(optionsDict(templateName))
{
    compileFiles_ = {templateName + "Template.C"};

    codedBase::setFilterVariable("typeName", codeName());

    forAll(substitutions_, i)
    {
        setFilterVariable(substitutions_[i]);
    }

    // Make verbose if debugging
    codedBase::setFilterVariable("verbose", Foam::name(bool(debug)));

    this->updateLibrary(dict_);
}


// ************************************************************************* //
