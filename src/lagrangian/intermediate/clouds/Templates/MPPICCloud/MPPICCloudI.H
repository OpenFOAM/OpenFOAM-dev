/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013 OpenFOAM Foundation
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

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
inline const Foam::MPPICCloud<CloudType>&
Foam::MPPICCloud<CloudType>::cloudCopy() const
{
    return cloudCopyPtr_();
}


template<class CloudType>
inline const Foam::PackingModel<Foam::MPPICCloud<CloudType> >&
Foam::MPPICCloud<CloudType>::packingModel() const
{
    return packingModel_();
}


template<class CloudType>
inline Foam::PackingModel<Foam::MPPICCloud<CloudType> >&
Foam::MPPICCloud<CloudType>::packingModel()
{
    return packingModel_();
}


template<class CloudType>
inline const Foam::DampingModel<Foam::MPPICCloud<CloudType> >&
Foam::MPPICCloud<CloudType>::dampingModel() const
{
    return dampingModel_();
}


template<class CloudType>
inline Foam::DampingModel<Foam::MPPICCloud<CloudType> >&
Foam::MPPICCloud<CloudType>::dampingModel()
{
    return dampingModel_();
}


template<class CloudType>
inline const Foam::IsotropyModel<Foam::MPPICCloud<CloudType> >&
Foam::MPPICCloud<CloudType>::isotropyModel() const
{
    return isotropyModel_();
}


template<class CloudType>
inline Foam::IsotropyModel<Foam::MPPICCloud<CloudType> >&
Foam::MPPICCloud<CloudType>::isotropyModel()
{
    return isotropyModel_();
}


// ************************************************************************* //
