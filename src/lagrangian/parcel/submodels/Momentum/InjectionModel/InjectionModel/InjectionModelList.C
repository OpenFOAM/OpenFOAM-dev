/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2012-2023 OpenFOAM Foundation
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

#include "InjectionModel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::InjectionModelList<CloudType>::InjectionModelList(CloudType& owner)
:
    PtrList<InjectionModel<CloudType>>()
{}


template<class CloudType>
Foam::InjectionModelList<CloudType>::InjectionModelList
(
    const dictionary& dict,
    CloudType& owner
)
:
    PtrList<InjectionModel<CloudType>>()
{
    wordList modelNames(dict.toc());

    Info<< "Constructing particle injection models" << endl;

    if (modelNames.size() > 0)
    {
        this->setSize(modelNames.size());

        label i = 0;
        forAllConstIter(IDLList<entry>, dict, iter)
        {
            const word& model = iter().keyword();
            Info<< "Creating injector: " << model << endl;
            const dictionary& props = iter().dict();

            this->set
            (
                i++,
                InjectionModel<CloudType>::New
                (
                    props,
                    model,
                    props.lookup("type"),
                    owner
                )
            );
        }
    }
    else
    {
        this->setSize(1);

        this->set
        (
            0,
            InjectionModel<CloudType>::New
            (
                dict,
                "none",
                "none",
                owner
            )
        );
    }
}


template<class CloudType>
Foam::InjectionModelList<CloudType>::InjectionModelList
(
    const InjectionModelList<CloudType>& iml
)
:
    PtrList<InjectionModel<CloudType>>(iml)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::InjectionModelList<CloudType>::~InjectionModelList()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
Foam::scalar Foam::InjectionModelList<CloudType>::timeStart() const
{
    scalar minTime = great;
    forAll(*this, i)
    {
        minTime = min(minTime, this->operator[](i).timeStart());
    }

    return minTime;
}


template<class CloudType>
Foam::scalar Foam::InjectionModelList<CloudType>::timeEnd() const
{
    scalar maxTime = -great;
    forAll(*this, i)
    {
        maxTime = max(maxTime, this->operator[](i).timeEnd());
    }

    return maxTime;
}


template<class CloudType>
void Foam::InjectionModelList<CloudType>::topoChange()
{
    forAll(*this, i)
    {
        this->operator[](i).topoChange();
    }
}


template<class CloudType>
template<class TrackCloudType>
void Foam::InjectionModelList<CloudType>::inject
(
    TrackCloudType& cloud,
    typename CloudType::parcelType::trackingData& td
)
{
    forAll(*this, i)
    {
        this->operator[](i).inject(cloud, td);
    }
}


template<class CloudType>
template<class TrackCloudType>
void Foam::InjectionModelList<CloudType>::injectSteadyState
(
    TrackCloudType& cloud,
    typename CloudType::parcelType::trackingData& td
)
{
    forAll(*this, i)
    {
        this->operator[](i).injectSteadyState(cloud, td);
    }
}


template<class CloudType>
void Foam::InjectionModelList<CloudType>::info(Ostream& os)
{
    forAll(*this, i)
    {
        this->operator[](i).info(os);
    }
}


// ************************************************************************* //
