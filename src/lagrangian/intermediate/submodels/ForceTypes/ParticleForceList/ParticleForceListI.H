/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class CloudType>
inline const CloudType& Foam::ParticleForceList<CloudType>::owner() const
{
    return owner_;
}


template<class CloudType>
inline CloudType& Foam::ParticleForceList<CloudType>::owner()
{
    return owner_;
}


template<class CloudType>
inline const Foam::fvMesh& Foam::ParticleForceList<CloudType>::mesh() const
{
    return mesh_;
}


template<class CloudType>
inline const Foam::dictionary& Foam::ParticleForceList<CloudType>::dict() const
{
    return dict_;
}


template<class CloudType>
inline void Foam::ParticleForceList<CloudType>::setCalcCoupled(bool flag)
{
    calcCoupled_ = flag;
}


template<class CloudType>
inline void Foam::ParticleForceList<CloudType>::setCalcNonCoupled(bool flag)
{
    calcNonCoupled_ = flag;
}


// ************************************************************************* //
