/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2015-2022 OpenFOAM Foundation
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

#include "diffusiveMassTransferModel.H"

// * * * * * * * * * * * * * * * * Selector  * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::diffusiveMassTransferModel>
Foam::diffusiveMassTransferModel::New
(
    const dictionary& dict,
    const phaseInterface& interface
)
{
    const dictionary& modelDict =
        interface.fluid().modelSubDict<diffusiveMassTransferModel>(dict);

    const word diffusiveMassTransferModelType(modelDict.lookup("type"));

    Info<< "Selecting diffusiveMassTransferModel for "
        << interface.name() << ": " << diffusiveMassTransferModelType << endl;

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(diffusiveMassTransferModelType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown diffusiveMassTransferModelType type "
            << diffusiveMassTransferModelType << endl << endl
            << "Valid diffusiveMassTransferModel types are : " << endl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return cstrIter()(modelDict, interface);
}


Foam::autoPtr<Foam::blendedDiffusiveMassTransferModel>
Foam::blendedDiffusiveMassTransferModel::New
(
    const dictionary& dict,
    const phaseInterface& interface
)
{
    return autoPtr<blendedDiffusiveMassTransferModel>
    (
        new blendedDiffusiveMassTransferModel(dict, interface)
    );
}


Foam::autoPtr<Foam::sidedBlendedDiffusiveMassTransferModel>
Foam::sidedBlendedDiffusiveMassTransferModel::New
(
    const dictionary& dict,
    const phaseInterface& interface
)
{
    return autoPtr<sidedBlendedDiffusiveMassTransferModel>
    (
        new sidedBlendedDiffusiveMassTransferModel(dict, interface)
    );
}


// ************************************************************************* //
