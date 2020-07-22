/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2013-2020 OpenFOAM Foundation
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

#include "IsotropyModel.H"

#include "TimeScaleModel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::IsotropyModel<CloudType>::IsotropyModel(CloudType& owner)
:
    CloudSubModelBase<CloudType>(owner),
    timeScaleModel_(nullptr)
{}


template<class CloudType>
Foam::IsotropyModel<CloudType>::IsotropyModel
(
    const dictionary& dict,
    CloudType& owner,
    const word& type
)
:
    CloudSubModelBase<CloudType>(owner, dict, typeName, type),
    timeScaleModel_
    (
        TimeScaleModel::New
        (
            this->coeffDict().subDict(TimeScaleModel::typeName)
        )
    )
{}


template<class CloudType>
Foam::IsotropyModel<CloudType>::IsotropyModel
(
    const IsotropyModel<CloudType>& cm
)
:
    CloudSubModelBase<CloudType>(cm),
    timeScaleModel_(cm.timeScaleModel_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::IsotropyModel<CloudType>::~IsotropyModel()
{}


// * * * * * * * * * * * * * * * *  Selector * * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::autoPtr<Foam::IsotropyModel<CloudType>>
Foam::IsotropyModel<CloudType>::New
(
    const dictionary& dict,
    CloudType& owner
)
{
    word modelType(dict.lookup(typeName));

    Info<< "Selecting isotropy model " << modelType << endl;

    typename dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(modelType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown isotropy model type " << modelType
            << ", constructor not in hash table" << nl << nl
            << "    Valid isotropy model types are:" << nl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return
        autoPtr<IsotropyModel<CloudType>>
        (
            cstrIter()(dict, owner)
        );
}


// ************************************************************************* //
