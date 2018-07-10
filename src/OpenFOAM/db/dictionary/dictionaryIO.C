/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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

#include "dictionary.H"
#include "inputModeEntry.H"
#include "regExp.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::dictionary::dictionary
(
    const fileName& name,
    const dictionary& parentDict,
    Istream& is
)
:
    dictionaryName(parentDict.name() + '.' + name),
    parent_(parentDict)
{
    read(is);
}


Foam::dictionary::dictionary(Istream& is)
:
    dictionaryName(is.name()),
    parent_(dictionary::null)
{
    // Reset input mode as this is a "top-level" dictionary
    functionEntries::inputModeEntry::clear();

    read(is);
}


Foam::dictionary::dictionary(Istream& is, const bool keepHeader)
:
    dictionaryName(is.name()),
    parent_(dictionary::null)
{
    // Reset input mode as this is a "top-level" dictionary
    functionEntries::inputModeEntry::clear();

    read(is, keepHeader);
}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::dictionary> Foam::dictionary::New(Istream& is)
{
    return autoPtr<dictionary>(new dictionary(is));
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::dictionary::read(Istream& is, const bool keepHeader)
{
    // Check for empty dictionary
    if (is.eof())
    {
        return true;
    }

    if (!is.good())
    {
        FatalIOErrorInFunction(is)
            << "Istream not OK for reading dictionary "
            << exit(FatalIOError);

        return false;
    }

    token currToken(is);
    if (currToken != token::BEGIN_BLOCK)
    {
        is.putBack(currToken);
    }

    while (!is.eof() && entry::New(*this, is))
    {}

    // normally remove the FoamFile header entry if it exists
    if (!keepHeader)
    {
        remove("FoamFile");
    }

    if (is.bad())
    {
        InfoInFunction
            << "Istream not OK after reading dictionary " << name()
            << endl;

        return false;
    }

    return true;
}


bool Foam::dictionary::read(Istream& is)
{
    return this->read(is, false);
}


bool Foam::dictionary::substituteKeyword(const word& keyword)
{
    word varName = keyword(1, keyword.size()-1);

    // lookup the variable name in the given dictionary
    const entry* ePtr = lookupEntryPtr(varName, true, true);

    // if defined insert its entries into this dictionary
    if (ePtr != nullptr)
    {
        const dictionary& addDict = ePtr->dict();

        forAllConstIter(IDLList<entry>, addDict, iter)
        {
            add(iter());
        }

        return true;
    }

    return false;
}


// * * * * * * * * * * * * * * Istream Operator  * * * * * * * * * * * * * * //

Foam::Istream& Foam::operator>>(Istream& is, dictionary& dict)
{
    // Reset input mode assuming this is a "top-level" dictionary
    functionEntries::inputModeEntry::clear();

    dict.clear();
    dict.name() = is.name();
    dict.read(is);

    return is;
}


// * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * * //

void Foam::dictionary::write(Ostream& os, bool subDict) const
{
    if (subDict)
    {
        os  << nl << indent << token::BEGIN_BLOCK << incrIndent << nl;
    }

    forAllConstIter(IDLList<entry>, *this, iter)
    {
        const entry& e = *iter;

        // Write entry
        os  << e;

        // Add extra new line between entries for "top-level" dictionaries
        if (!subDict && parent() == dictionary::null && e != *last())
        {
            os  << nl;
        }

        // Check stream before going to next entry.
        if (!os.good())
        {
            WarningInFunction
                << "Can't write entry " << iter().keyword()
                << " for dictionary " << name()
                << endl;
        }
    }

    if (subDict)
    {
        os  << decrIndent << indent << token::END_BLOCK << endl;
    }
}


Foam::Ostream& Foam::operator<<(Ostream& os, const dictionary& dict)
{
    dict.write(os, true);
    return os;
}


// ************************************************************************* //
