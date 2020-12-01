/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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

#include "CloudSubModelBase.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::CloudSubModelBase<CloudType>::CloudSubModelBase(CloudType& owner)
:
    subModelBase(owner.outputProperties()),
    owner_(owner)
{}


template<class CloudType>
Foam::CloudSubModelBase<CloudType>::CloudSubModelBase
(
    CloudType& owner,
    const dictionary& dict,
    const word& baseName,
    const word& modelType,
    const word& dictExt
)
:
    subModelBase
    (
        owner.outputProperties(),
        dict,
        baseName,
        modelType,
        dictExt
    ),
    owner_(owner)
{}


template<class CloudType>
Foam::CloudSubModelBase<CloudType>::CloudSubModelBase
(
    const word& modelName,
    CloudType& owner,
    const dictionary& dict,
    const word& baseName,
    const word& modelType
)
:
    subModelBase
    (
        modelName,
        owner.outputProperties(),
        dict,
        baseName,
        modelType
    ),
    owner_(owner)
{}


template<class CloudType>
Foam::CloudSubModelBase<CloudType>::CloudSubModelBase
(
    const CloudSubModelBase<CloudType>& smb
)
:
    subModelBase(smb),
    owner_(smb.owner_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::CloudSubModelBase<CloudType>::~CloudSubModelBase()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
const CloudType& Foam::CloudSubModelBase<CloudType>::owner() const
{
    return owner_;
}


template<class CloudType>
CloudType& Foam::CloudSubModelBase<CloudType>::owner()
{
    return owner_;
}


template<class CloudType>
bool Foam::CloudSubModelBase<CloudType>::writeTime() const
{
    return
        owner_.solution().transient()
     && owner_.db().time().writeTime();
}


template<class CloudType>
void Foam::CloudSubModelBase<CloudType>::write(Ostream& os) const
{
    writeEntry(os, "owner", owner_.name());
    subModelBase::write(os);
}


// ************************************************************************* //
