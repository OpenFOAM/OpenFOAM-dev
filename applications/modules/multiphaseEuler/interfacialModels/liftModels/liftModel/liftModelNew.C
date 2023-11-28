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

#include "liftModel.H"
#include "phaseSystem.H"

// * * * * * * * * * * * * * * * * Selector  * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::liftModel> Foam::liftModel::New
(
    const dictionary& dict,
    const phaseInterface& interface,
    const bool outer
)
{
    const dictionary& modelDict =
        outer ? interface.fluid().modelSubDict<liftModel>(dict) : dict;

    const word liftModelType(modelDict.lookup("type"));

    Info<< "Selecting liftModel for "
        << interface.name() << ": " << liftModelType << endl;

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(liftModelType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown liftModel type "
            << liftModelType << endl << endl
            << "Valid liftModel types are : " << endl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return cstrIter()(modelDict, interface);
}


Foam::autoPtr<Foam::blendedLiftModel> Foam::blendedLiftModel::New
(
    const dictionary& dict,
    const phaseInterface& interface
)
{
    return autoPtr<blendedLiftModel>(new blendedLiftModel(dict, interface));
}


// ************************************************************************* //
