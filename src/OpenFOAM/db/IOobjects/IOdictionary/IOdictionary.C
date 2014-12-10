/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
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

#include "IOdictionary.H"
#include "objectRegistry.H"
#include "Pstream.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
defineTypeNameAndDebug(IOdictionary, 0);

bool IOdictionary::writeDictionaries
(
    debug::infoSwitch("writeDictionaries", 0)
);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::IOdictionary::IOdictionary(const IOobject& io)
:
    regIOobject(io)
{
    // Temporary warning
    if (debug && io.readOpt() == IOobject::MUST_READ)
    {
        WarningIn("IOdictionary::IOdictionary(const IOobject&)")
            << "Dictionary " << name()
            << " constructed with IOobject::MUST_READ"
            " instead of IOobject::MUST_READ_IF_MODIFIED." << nl
            << "Use MUST_READ_IF_MODIFIED if you need automatic rereading."
            << endl;
    }

    // Everyone check or just master
    bool masterOnly =
        regIOobject::fileModificationChecking == timeStampMaster
     || regIOobject::fileModificationChecking == inotifyMaster;


    // Check if header is ok for READ_IF_PRESENT
    bool isHeaderOk = false;
    if (io.readOpt() == IOobject::READ_IF_PRESENT)
    {
        if (masterOnly)
        {
            if (Pstream::master())
            {
                isHeaderOk = headerOk();
            }
            Pstream::scatter(isHeaderOk);
        }
        else
        {
            isHeaderOk = headerOk();
        }
    }


    if
    (
        (
            io.readOpt() == IOobject::MUST_READ
         || io.readOpt() == IOobject::MUST_READ_IF_MODIFIED
        )
     || isHeaderOk
    )
    {
        readFile(masterOnly);
    }

    dictionary::name() = IOobject::objectPath();
}


Foam::IOdictionary::IOdictionary(const IOobject& io, const dictionary& dict)
:
    regIOobject(io)
{
    // Temporary warning
    if (debug && io.readOpt() == IOobject::MUST_READ)
    {
        WarningIn
        (
            "IOdictionary::IOdictionary(const IOobject& const dictionary&)"
        )   << "Dictionary " << name()
            << " constructed with IOobject::MUST_READ"
            " instead of IOobject::MUST_READ_IF_MODIFIED." << nl
            << "Use MUST_READ_IF_MODIFIED if you need automatic rereading."
            << endl;
    }

    // Everyone check or just master
    bool masterOnly =
        regIOobject::fileModificationChecking == timeStampMaster
     || regIOobject::fileModificationChecking == inotifyMaster;


    // Check if header is ok for READ_IF_PRESENT
    bool isHeaderOk = false;
    if (io.readOpt() == IOobject::READ_IF_PRESENT)
    {
        if (masterOnly)
        {
            if (Pstream::master())
            {
                isHeaderOk = headerOk();
            }
            Pstream::scatter(isHeaderOk);
        }
        else
        {
            isHeaderOk = headerOk();
        }
    }


    if
    (
        (
            io.readOpt() == IOobject::MUST_READ
         || io.readOpt() == IOobject::MUST_READ_IF_MODIFIED
        )
     || isHeaderOk
    )
    {
        readFile(masterOnly);
    }
    else
    {
        dictionary::operator=(dict);
    }

    dictionary::name() = IOobject::objectPath();
}


Foam::IOdictionary::IOdictionary(const IOobject& io, Istream& is)
:
    regIOobject(io)
{
    dictionary::name() = IOobject::objectPath();
    // Note that we do construct the dictionary null and read in afterwards
    // so that if there is some fancy massaging due to a functionEntry in
    // the dictionary at least the type information is already complete.
    is  >> *this;
}


// * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * * //

Foam::IOdictionary::~IOdictionary()
{}


// * * * * * * * * * * * * * * * Members Functions * * * * * * * * * * * * * //

const Foam::word& Foam::IOdictionary::name() const
{
    return regIOobject::name();
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

void Foam::IOdictionary::operator=(const IOdictionary& rhs)
{
    dictionary::operator=(rhs);
}


// ************************************************************************* //
