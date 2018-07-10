/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2013-2018 OpenFOAM Foundation
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
#include "Time.H"
#include "polyMesh.H"
#include "addToRunTimeSelectionTable.H"
#include "IOdictionary.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(writeDictionary, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        writeDictionary,
        dictionary
    );
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::functionObjects::writeDictionary::tryDirectory
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

    if (dictIO.typeHeaderOk<IOdictionary>(true))
    {
        IOdictionary dict(dictIO);

        if (dict.digest() != digests_[dictI])
        {
            if (firstDict)
            {
                Info<< type() << " " << name() << " write:" << nl << endl;

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

Foam::functionObjects::writeDictionary::writeDictionary
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    functionObject(name),
    obr_
    (
        runTime.lookupObject<objectRegistry>
        (
            dict.lookupOrDefault("region", polyMesh::defaultRegion)
        )
    ),
    dictNames_(),
    digests_()
{
    read(dict);
    execute();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::writeDictionary::~writeDictionary()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::writeDictionary::read(const dictionary& dict)
{
    wordList dictNames(dict.lookup("dictNames"));
    HashSet<word> uniqueNames(dictNames);
    dictNames_ = uniqueNames.toc();

    digests_.setSize(dictNames_.size(), SHA1Digest());

    Info<< type() << " " << name() << ": monitoring dictionaries:" << nl;
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

    return true;
}


bool Foam::functionObjects::writeDictionary::execute()
{
    return true;
}


bool Foam::functionObjects::writeDictionary::write()
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
                    Info<< type() << " " << name() << " write:" << nl << endl;

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

    return true;
}


// ************************************************************************* //
