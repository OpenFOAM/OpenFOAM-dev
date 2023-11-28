/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2014-2023 OpenFOAM Foundation
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

#include "wallLubricationModel.H"
#include "phaseSystem.H"

// * * * * * * * * * * * * * * * * Selector  * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::wallLubricationModel> Foam::wallLubricationModel::New
(
    const dictionary& dict,
    const phaseInterface& interface,
    const bool outer
)
{
    const dictionary& modelDict =
        outer
      ? interface.fluid().modelSubDict<wallLubricationModel>(dict)
      : dict;

    const word wallLubricationModelType(modelDict.lookup("type"));

    Info<< "Selecting wallLubricationModel for "
        << interface.name() << ": " << wallLubricationModelType << endl;

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(wallLubricationModelType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown wallLubricationModel type "
            << wallLubricationModelType << endl << endl
            << "Valid wallLubricationModel types are : " << endl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return cstrIter()(modelDict, interface);
}


Foam::autoPtr<Foam::blendedWallLubricationModel>
Foam::blendedWallLubricationModel::New
(
    const dictionary& dict,
    const phaseInterface& interface
)
{
    return
        autoPtr<blendedWallLubricationModel>
        (
            new blendedWallLubricationModel(dict, interface)
        );
}


// ************************************************************************* //
