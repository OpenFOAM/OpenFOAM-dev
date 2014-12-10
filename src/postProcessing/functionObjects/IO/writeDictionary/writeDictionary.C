/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013-2014 OpenFOAM Foundation
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

#include "writeDictionary.H"
#include "dictionary.H"
#include "Time.H"
#include "HashSet.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
defineTypeNameAndDebug(writeDictionary, 0);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::writeDictionary::tryDirectory
(
    const label dictI,
    const word& location,
    bool& firstDict
)
{
    IOobject dictIO
    (
        dictNames_[dictI],
        location,
        obr_,
        IOobject::MUST_READ,
        IOobject::NO_WRITE,
        false
    );

    if (dictIO.headerOk())
    {
        IOdictionary dict(dictIO);

        if (dict.digest() != digests_[dictI])
        {
            if (firstDict)
            {
                Info<< type() << " " << name_ << " output:" << nl << endl;

                IOobject::writeDivider(Info);
                Info<< endl;
                firstDict = false;
            }

            Info<< dict.dictName() << dict << nl;

            IOobject::writeDivider(Info);

            digests_[dictI] = dict.digest();
        }

        return true;
    }

    return false;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::writeDictionary::writeDictionary
(
    const word& name,
    const objectRegistry& obr,
    const dictionary& dict,
    const bool loadFromFiles
)
:
    name_(name),
    obr_(obr),
    dictNames_(),
    digests_()
{
    read(dict);

    execute();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::writeDictionary::~writeDictionary()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::writeDictionary::read(const dictionary& dict)
{
    wordList dictNames(dict.lookup("dictNames"));
    HashSet<word> uniqueNames(dictNames);
    dictNames_ = uniqueNames.toc();

    digests_.setSize(dictNames_.size(), SHA1Digest());

    Info<< type() << " " << name_ << ": monitoring dictionaries:" << nl;
    if (dictNames_.size())
    {
        forAll(dictNames_, i)
        {
            Info<< "    " << dictNames_[i] << endl;
        }
    }
    else
    {
        Info<< "    none" << nl;
    }
    Info<< endl;
}


void Foam::writeDictionary::execute()
{
    bool firstDict = true;
    forAll(dictNames_, i)
    {
        if (obr_.foundObject<dictionary>(dictNames_[i]))
        {
            const dictionary& dict =
                obr_.lookupObject<dictionary>(dictNames_[i]);

            if (dict.digest() != digests_[i])
            {
                if (firstDict)
                {
                    Info<< type() << " " << name_ << " output:" << nl << endl;

                    IOobject::writeDivider(Info);
                    Info<< endl;
                    firstDict = false;
                }

                digests_[i] = dict.digest();

                Info<< dict.dictName() << dict << nl;
                IOobject::writeDivider(Info);
                Info<< endl;
            }
        }
        else
        {
            bool processed = tryDirectory(i, obr_.time().timeName(), firstDict);

            if (!processed)
            {
                processed = tryDirectory(i, obr_.time().constant(), firstDict);
            }

            if (!processed)
            {
                processed = tryDirectory(i, obr_.time().system(), firstDict);
            }

            if (!processed)
            {
                Info<< "    Unable to locate dictionary " << dictNames_[i]
                    << nl << endl;
            }
            else
            {
                Info<< endl;
            }
        }
    }
}


void Foam::writeDictionary::end()
{
    execute();
}


void Foam::writeDictionary::timeSet()
{
    // do nothing
}


void Foam::writeDictionary::write()
{
    // do nothing
}


// ************************************************************************* //
