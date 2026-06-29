/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2026 OpenFOAM Foundation
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

#include "solidProperties.H"
#include "Switch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::solidProperties> Foam::solidProperties::New
(
    const word& name
)
{
    if (debug)
    {
        InfoInFunction << "Constructing solidProperties" << endl;
    }

    ConstructorTable::iterator cstrIter = ConstructorTablePtr_->find(name);

    if (cstrIter == ConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown solidProperties type "
            << name << nl << nl
            << "Valid solidProperties types are:" << nl
            << ConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<solidProperties>(cstrIter()());
}


Foam::autoPtr<Foam::solidProperties> Foam::solidProperties::New
(
    const dictionary& dict
)
{
    if (debug)
    {
        InfoInFunction << "Constructing solid" << endl;
    }

    // If the type is not specified use the dict name as the liquid type name
    const word& solidPropertiesTypeName =
        dict.found("type") ? dict.lookup<word>("type") : dict.dictName();

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(solidPropertiesTypeName);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalIOErrorInFunction(dict)
            << "Unknown solidProperties type "
            << solidPropertiesTypeName << nl << nl
            << "Valid solidProperties types are:" << nl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalIOError);
    }

    return autoPtr<solidProperties>(cstrIter()(dict));
}


// ************************************************************************* //
