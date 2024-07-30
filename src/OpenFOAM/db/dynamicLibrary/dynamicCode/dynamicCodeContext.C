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
    const dictionary& contextDict,
    const dictionary& codeDict,
    const wordList& codeKeys,
    const wordList& codeDictVars
)
:
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
        codePtrs[i] = codeDict.lookupEntryPtr(key, false, false);
        if (codePtrs[i])
        {
            string s(stringOps::trim(verbatimString(codePtrs[i]->stream())));
            stringOps::inplaceExpandCodeString
            (
                s,
                contextDict, // Lookup variables from the context dictionary
                codeDictVars[i]
            );
            code_.insert(key, s);
        }
        else
        {
            code_.insert(key, "");
        }
    }

    // Options
    const entry* optionsPtr =
        codeDict.lookupEntryPtr("codeOptions", false, false);
    if (optionsPtr)
    {
        options_ = stringOps::trim(verbatimString(optionsPtr->stream()));
        stringOps::inplaceExpandCodeString(options_, contextDict, word::null);
    }

    // Libs
    const entry* libsPtr = codeDict.lookupEntryPtr("codeLibs", false, false);
    if (libsPtr)
    {
        libs_ = stringOps::trim(verbatimString(libsPtr->stream()));
        stringOps::inplaceExpandCodeString(libs_, contextDict, word::null);
    }

    // Calculate SHA1 digest from all entries
    OSHA1stream os;
    forAllConstIter(HashTable<string>, code_, iter)
    {
        os << iter();
    }
    os << options_ << libs_;
    sha1_ = os.digest();

    // Add line directives after calculating SHA1
    forAll(codeKeys, i)
    {
        if (codePtrs[i])
        {
            const word& key = codeKeys[i];
            addLineDirective
            (
                code_[key],
                codePtrs[i]->startLineNumber(),
                codeDict.name()
            );
        }
    }
}


Foam::dynamicCodeContext::dynamicCodeContext
(
    const dictionary& contextDict,
    const wordList& codeKeys,
    const wordList& codeDictVars
)
:
    dynamicCodeContext(contextDict, contextDict, codeKeys, codeDictVars)
{}


// ************************************************************************* //
