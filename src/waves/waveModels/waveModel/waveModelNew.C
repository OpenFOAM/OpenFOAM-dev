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

#include "waveModel.H"

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::waveModel> Foam::waveModel::New
(
    const objectRegistry& db,
    const dictionary& dict
)
{
    return waveModel::New(dict.lookup("type"), db, dict);
}


Foam::autoPtr<Foam::waveModel> Foam::waveModel::New
(
    const word& type,
    const objectRegistry& db,
    const dictionary& dict
)
{
    if (debug)
    {
        Info<< "Selecting " << waveModel::typeName << " " << type << endl;
    }

    objectRegistryConstructorTable::iterator cstrIter =
        objectRegistryConstructorTablePtr_->find(type);

    if (cstrIter == objectRegistryConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown " << waveModel::typeName << " " << type << nl << nl
            << "Valid model types are:" << nl
            << objectRegistryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return cstrIter()(db, dict);
}


// ************************************************************************* //
