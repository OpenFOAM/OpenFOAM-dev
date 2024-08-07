/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2024 OpenFOAM Foundation
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

#include "heatTransferModel.H"
#include "phaseSystem.H"

// * * * * * * * * * * * * * * * * Selector  * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::heatTransferModel> Foam::heatTransferModel::New
(
    const dictionary& dict,
    const phaseInterface& interface,
    const bool outer,
    const bool registerObject
)
{
    const dictionary& modelDict =
        outer ? interface.fluid().modelSubDict<heatTransferModel>(dict) : dict;

    const word heatTransferModelType(modelDict.lookup("type"));

    Info<< "Selecting heatTransferModel for "
        << interface.name() << ": " << heatTransferModelType << endl;

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(heatTransferModelType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalIOErrorInFunction(dict)
            << "Unknown heatTransferModel type "
            << heatTransferModelType << endl << endl
            << "Valid heatTransferModel types are : " << endl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalIOError);
    }

    return cstrIter()(modelDict, interface, registerObject);
}


Foam::autoPtr<Foam::blendedHeatTransferModel>
Foam::blendedHeatTransferModel::New
(
    const dictionary& dict,
    const phaseInterface& interface
)
{
    return autoPtr<blendedHeatTransferModel>
    (
        new blendedHeatTransferModel(dict, interface)
    );
}


Foam::autoPtr<Foam::sidedBlendedHeatTransferModel>
Foam::sidedBlendedHeatTransferModel::New
(
    const dictionary& dict,
    const phaseInterface& interface
)
{
    return autoPtr<sidedBlendedHeatTransferModel>
    (
        new sidedBlendedHeatTransferModel(dict, interface)
    );
}


// ************************************************************************* //
