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


void Foam::dynamicCodeContext::read
(
    const dictionary& contextDict,
    const dictionary& codeDict
)
{
    // Expand all dictionary entries. Note that this removes any leading or
    // trailing whitespace, which is necessary for compilation options, and
    // doesn't hurt for everything else
    List<const entry*> codePtrs(codeKeys_.size(), nullptr);
    code_.clear();
    forAll(codeKeys_, i)
    {
        const word& key = codeKeys_[i];
        codePtrs[i] = codeDict.lookupEntryPtr(key, false, false);
        if (codePtrs[i])
        {
            codeStrings_[i] = verbatimString(codePtrs[i]->stream());
            stringOps::inplaceExpandCodeString
            (
                codeStrings_[i],
                contextDict, // Lookup variables from the context dictionary
                codeDictVars_[i]
            );
            code_.insert(key, codeStrings_[i]);
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
        optionsString_ = verbatimString(optionsPtr->stream());
        stringOps::inplaceExpandCodeString
        (
            optionsString_,
            contextDict,
            word::null
        );
        options_ = stringOps::trim(optionsString_);
    }

    // Libs
    const entry* libsPtr = codeDict.lookupEntryPtr("codeLibs", false, false);
    if (libsPtr)
    {
        libsString_ = verbatimString(libsPtr->stream());
        stringOps::inplaceExpandCodeString
        (
            libsString_,
            contextDict,
            word::null
        );
        libs_ = stringOps::trim(libsString_);
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
    forAll(codeKeys_, i)
    {
        if (codePtrs[i])
        {
            const word& key = codeKeys_[i];
            addLineDirective
            (
                code_[key],
                codePtrs[i]->startLineNumber(),
                codeDict.name()
            );
        }
    }
}


void Foam::dynamicCodeContext::read(const dictionary& contextDict)
{
    read(contextDict, contextDict);
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
    codeKeys_(codeKeys),
    codeDictVars_(codeDictVars),
    codeStrings_(codeKeys.size())
{
    read(contextDict, codeDict);
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


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::dynamicCodeContext::write(Ostream& os) const
{
    forAll(codeStrings_, i)
    {
        if (codeStrings_[i] != verbatimString::null)
        {
            writeEntry(os, codeKeys_[i], codeStrings_[i]);
        }
    }

    if (optionsString_ != verbatimString::null)
    {
        writeEntry(os, "codeOptions", optionsString_);
    }

    if (libsString_ != verbatimString::null)
    {
        writeEntry(os, "codeLibs", libsString_);
    }
}


// ************************************************************************* //
