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

#include "PackingModel.H"
#include "AveragingMethod.H"
#include "ParticleStressModel.H"
#include "CorrectionLimitingMethod.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::PackingModel<CloudType>::PackingModel(CloudType& owner)
:
    CloudSubModelBase<CloudType>(owner),
    particleStressModel_(nullptr)
{}


template<class CloudType>
Foam::PackingModel<CloudType>::PackingModel
(
    const dictionary& dict,
    CloudType& owner,
    const word& type
)
:
    CloudSubModelBase<CloudType>(owner, dict, typeName, type),
    particleStressModel_
    (
        ParticleStressModel::New
        (
            this->coeffDict().subDict(ParticleStressModel::typeName)
        )
    )
{}


template<class CloudType>
Foam::PackingModel<CloudType>::PackingModel(const PackingModel<CloudType>& cm)
:
    CloudSubModelBase<CloudType>(cm),
    particleStressModel_(cm.particleStressModel_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::PackingModel<CloudType>::~PackingModel()
{}


// * * * * * * * * * * * * * * * * Selector  * * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::autoPtr<Foam::PackingModel<CloudType>>
Foam::PackingModel<CloudType>::New
(
    const dictionary& dict,
    CloudType& owner
)
{
    word modelType(dict.lookup(typeName));

    Info<< "Selecting packing model " << modelType << endl;

    typename dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(modelType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown packing model type " << modelType
            << ", constructor not in hash table" << nl << nl
            << "    Valid packing model types are:" << nl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<PackingModel<CloudType>>(cstrIter()(dict, owner));
}


// ************************************************************************* //
