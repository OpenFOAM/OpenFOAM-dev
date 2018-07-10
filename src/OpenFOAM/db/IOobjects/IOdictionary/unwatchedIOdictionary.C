/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2015-2018 OpenFOAM Foundation
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

#include "unwatchedIOdictionary.H"
#include "objectRegistry.H"
#include "Pstream.H"
#include "Time.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::unwatchedIOdictionary::unwatchedIOdictionary(const IOobject& io)
:
    baseIOdictionary(io)
{
    readHeaderOk(IOstream::ASCII, typeName);

    // For if MUST_READ_IF_MODIFIED
    addWatch();
}


Foam::unwatchedIOdictionary::unwatchedIOdictionary
(
    const IOobject& io,
    const dictionary& dict
)
:
    baseIOdictionary(io, dict)
{
    if (!readHeaderOk(IOstream::ASCII, typeName))
    {
        dictionary::operator=(dict);
    }

    // For if MUST_READ_IF_MODIFIED
    addWatch();
}


Foam::unwatchedIOdictionary::unwatchedIOdictionary
(
    const IOobject& io,
    Istream& is
)
:
    baseIOdictionary(io, is)
{
    // Note that we do construct the dictionary null and read in
    // afterwards
    // so that if there is some fancy massaging due to a
    // functionEntry in
    // the dictionary at least the type information is already complete.
    is  >> *this;

    // For if MUST_READ_IF_MODIFIED
    addWatch();
}


// * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * * //

Foam::unwatchedIOdictionary::~unwatchedIOdictionary()
{}


// * * * * * * * * * * * * * * * Members Functions * * * * * * * * * * * * * //

Foam::label Foam::unwatchedIOdictionary::addWatch(const fileName& f)
{
    label index = -1;

    if (readOpt() == MUST_READ_IF_MODIFIED)
    {
        index = findIndex(files_, f);

        if (index == -1)
        {
            index = files_.size();
            files_.append(f);
        }
    }
    return index;
}


void Foam::unwatchedIOdictionary::addWatch()
{
    if (readOpt() == MUST_READ_IF_MODIFIED)
    {
        fileName f = filePath();
        if (!f.size())
        {
            // We don't have this file but would like to re-read it.
            // Possibly if master-only reading mode.
            f = objectPath();
        }

        if (findIndex(files_, f) != -1)
        {
            FatalErrorIn("regIOobject::addWatch()")
                << "Object " << objectPath() << " of type " << type()
                << " already watched" << abort(FatalError);
        }

        // If master-only reading only the master will have all dependencies
        // so scatter these to slaves
        bool masterOnly =
            global()
         && (
                regIOobject::fileModificationChecking == timeStampMaster
             || regIOobject::fileModificationChecking == inotifyMaster
            );

        if (masterOnly && Pstream::parRun())
        {
            Pstream::scatter(files_);
        }

        addWatch(f);
    }
}


// ************************************************************************* //
