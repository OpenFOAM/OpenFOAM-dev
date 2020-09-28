/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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

#include "coordinateSystems.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::autoPtr<Foam::coordinateSystem> Foam::coordinateSystem::New
(
    const objectRegistry& obr,
    const dictionary& dict
)
{
    const dictionary& coordDict = dict.subDict(typeName_());

    const entry* entryPtr = coordDict.lookupEntryPtr(typeName_(), false, false);

    // non-dictionary entry is a lookup into global coordinateSystems
    if (entryPtr && !entryPtr->isDict())
    {
        keyType key(entryPtr->stream());

        const coordinateSystems& lst = coordinateSystems::New(obr);
        const label index = lst.findIndex(key);

        if (debug)
        {
            InfoInFunction
                << "Using global coordinate system: "
                << key << "=" << index << endl;
        }

        if (index < 0)
        {
            FatalErrorInFunction
                << "could not find coordinate system: " << key << nl
                << "available coordinate systems: " << lst.toc() << nl << nl
                << exit(FatalError);
        }

        return lst[index].clone();
    }
    else
    {
        const word coordType = coordDict.lookup("type");

        dictionaryConstructorTable::iterator cstrIter =
            dictionaryConstructorTablePtr_->find(coordType);

        if (cstrIter == dictionaryConstructorTablePtr_->end())
        {
            FatalIOErrorInFunction
            (
                dict
            )   << "Unknown coordinateSystem type "
                << coordType << nl << nl
                << "Valid coordinateSystem types are :" << nl
                << dictionaryConstructorTablePtr_->sortedToc()
                << exit(FatalIOError);
        }

        return autoPtr<coordinateSystem>(cstrIter()(coordType, coordDict));
    }
}


Foam::autoPtr<Foam::coordinateSystem> Foam::coordinateSystem::New
(
    const word& name,
    const dictionary& dict
)
{
    const word coordType = dict.lookup("type");

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(coordType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalIOErrorInFunction
        (
            dict
        )   << "Unknown coordinateSystem type "
            << coordType << nl << nl
            << "Valid coordinateSystem types are :" << nl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalIOError);
    }

    return autoPtr<coordinateSystem>(cstrIter()(name, dict));
}


Foam::autoPtr<Foam::coordinateSystem> Foam::coordinateSystem::New
(
    Istream& is
)
{
    const word name(is);
    const dictionary dict(is);

    return autoPtr<coordinateSystem>(coordinateSystem::New(name, dict));
}


// ************************************************************************* //
