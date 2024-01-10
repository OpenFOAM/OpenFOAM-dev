/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2024 OpenFOAM Foundation
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

#include "dynamicCodeContext.H"
#include "stringOps.H"
#include "OSHA1stream.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::dynamicCodeContext::addLineDirective
(
    string& code,
    const label lineNum,
    const fileName& name
)
{
    code = "#line " + Foam::name(lineNum + 1) + " \"" + name + "\"\n" + code;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::dynamicCodeContext::dynamicCodeContext
(
    const dictionary& dict,
    const wordList& codeKeys,
    const wordList& codeDictVars
)
:
    dict_(dict),
    code_(),
    options_(),
    libs_()
{
    // Expand all dictionary entries. Note that this removes any leading or
    // trailing whitespace, which is necessary for compilation options, and
    // doesn't hurt for everything else
    List<const entry*> codePtrs(codeKeys.size(), nullptr);
    forAll(codeKeys, i)
    {
        const word& key = codeKeys[i];
        codePtrs[i] = dict.lookupEntryPtr(key, false, false);
        if (codePtrs[i])
        {
            string s(stringOps::trim(verbatimString(codePtrs[i]->stream())));
            stringOps::inplaceExpandCodeString(s, dict, codeDictVars[i]);
            code_.insert(key, s);
        }
        else
        {
            code_.insert(key, "");
        }
    }

    // Options
    const entry* optionsPtr = dict.lookupEntryPtr("codeOptions", false, false);
    if (optionsPtr)
    {
        options_ = stringOps::trim(verbatimString(optionsPtr->stream()));
        stringOps::inplaceExpandCodeString(options_, dict, word::null);
    }

    // Libs
    const entry* libsPtr = dict.lookupEntryPtr("codeLibs", false, false);
    if (libsPtr)
    {
        libs_ = stringOps::trim(verbatimString(libsPtr->stream()));
        stringOps::inplaceExpandCodeString(libs_, dict, word::null);
    }

    // Calculate SHA1 digest from all entries
    OSHA1stream os;
    forAllConstIter(HashTable<string>, code_, iter)
    {
        os << iter();
    }
    os << options_ << libs_;
    sha1_ = os.digest();

    // Add line directive after calculating SHA1 since this includes
    // "processor..." in the path which differs between cores
    forAll(codeKeys, i)
    {
        if (codePtrs[i])
        {
            const word& key = codeKeys[i];
            addLineDirective
            (
                code_[key],
                codePtrs[i]->startLineNumber(),
                dict.name()
            );
        }
    }
}


// ************************************************************************* //
