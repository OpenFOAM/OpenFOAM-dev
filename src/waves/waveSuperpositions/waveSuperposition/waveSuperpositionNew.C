/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2017-2020 OpenFOAM Foundation
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

#include "waveSuperposition.H"
#include "Time.H"

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

const Foam::waveSuperposition& Foam::waveSuperposition::New
(
    const objectRegistry& db
)
{
    if (db.foundObject<waveSuperposition>(dictName))
    {
        return db.lookupObject<waveSuperposition>(dictName);
    }

    const IOdictionary dict
    (
        IOobject
        (
            dictName,
            db.time().constant(),
            db,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        )
    );

    const word type =
        dict.lookupOrDefault<word>("type", waveSuperposition::typeName);

    objectRegistryConstructorTable::iterator cstrIter =
        objectRegistryConstructorTablePtr_->find(type);

    if (cstrIter == objectRegistryConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown " << waveSuperposition::typeName << " " << type
            << nl << nl << "Valid types are:" << nl
            << objectRegistryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    waveSuperposition* ptr = cstrIter()(db).ptr();

    ptr->store();

    return *ptr;
}


// ************************************************************************* //
