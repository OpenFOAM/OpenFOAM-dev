/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2025 OpenFOAM Foundation
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

#include "codeBlockCalcEntry.H"
#include "codeBlockEntry.H"
#include "OSspecific.H"
#include "codeStream.H"
#include "addToRunTimeSelectionTable.H"
#include "addToMemberFunctionSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionEntries
{
    defineFunctionTypeNameAndDebug(codeBlockCalcEntry, 0);

    addToRunTimeSelectionTable(functionEntry, codeBlockCalcEntry, dictionary);

    addToMemberFunctionSelectionTable
    (
        functionEntry,
        codeBlockCalcEntry,
        execute,
        primitiveEntryIstream
    );
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::OTstream Foam::functionEntries::codeBlockCalcEntry::resultStream
(
    const dictionary& dict,
    Istream& is
)
{
    if (debug)
    {
        Info<< "Expanding " << typeName << " at line " << is.lineNumber()
            << " in file " <<  dict.name() << endl;
    }

    // Construct the name of the function
    // corresponding to this codeBlockCalcEntry
    const word functionName
    (
        codeBlockPtr_->codeBlockName_ + "_" + Foam::name(index_)
    );

    // Find the function handle in the library
    const codeStream::streamingFunctionType function =
        reinterpret_cast<codeStream::streamingFunctionType>
        (
            dlSym(codeBlockPtr_->lib_, functionName)
        );

    if (!function)
    {
        FatalIOErrorInFunction(dict)
            << "Failed looking up symbol " << functionName
            << " in library " << codeBlockPtr_->lib_ << exit(FatalIOError);
    }

    // Use function to write stream
    OTstream ots(is.name(), is.format());
    (*function)(ots, dict);

    // Return the OTstream containing the results of the calculation
    return ots;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionEntries::codeBlockCalcEntry::codeBlockCalcEntry
(
    const label lineNumber,
    const dictionary& parentDict,
    Istream& is
)
:
    functionEntry(typeName, lineNumber, parentDict),
    codeBlockPtr_(nullptr),
    index_(-1)
{
    token t(is);

    if (t.isUnsignedInteger64())
    {
        codeBlockPtr_ =
            reinterpret_cast<const codeBlockEntry*>(t.unsignedInteger64Token());

        is.read(t);

        if (t.isLabel())
        {
            index_ = t.labelToken();
        }
        else
        {
            FatalIOErrorInFunction(is)
                << "Wrong argument type for " << typeName << nl
                << "    Expected a label but found token " << t
                << exit(FatalIOError);
        }
    }
    else
    {
        FatalIOErrorInFunction(is)
            << "Wrong argument type for " << typeName << nl
            << "    Expected a unsignedInteger64 but found token " << t
            << exit(FatalIOError);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionEntries::codeBlockCalcEntry::execute
(
    dictionary& contextDict,
    Istream& is
)
{
    return insert(contextDict, resultStream(contextDict, is));
}


bool Foam::functionEntries::codeBlockCalcEntry::execute
(
    const dictionary& contextDict,
    primitiveEntry& contextEntry,
    Istream& is
)
{
    return insert
    (
        contextDict,
        contextEntry,
        codeBlockCalcEntry(is.lineNumber(), contextDict, is)
       .resultStream(contextDict, is)
    );
}


void Foam::functionEntries::codeBlockCalcEntry::write(Ostream& os) const
{}


// ************************************************************************* //
