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

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ParticleType>
const Foam::polyMesh& Foam::InteractionLists<ParticleType>::mesh() const
{
    return mesh_;
}


template<class ParticleType>
const Foam::mapDistribute&
Foam::InteractionLists<ParticleType>::cellMap() const
{
    return cellMapPtr_();
}


template<class ParticleType>
const Foam::mapDistribute&
Foam::InteractionLists<ParticleType>::wallFaceMap() const
{
    return wallFaceMapPtr_();
}


template<class ParticleType>
const Foam::labelListList& Foam::InteractionLists<ParticleType>::dil() const
{
    return dil_;
}


template<class ParticleType>
const Foam::labelListList&
Foam::InteractionLists<ParticleType>::dwfil() const
{
    return dwfil_;
}


template<class ParticleType>
const Foam::labelListList& Foam::InteractionLists<ParticleType>::ril() const
{
    return ril_;
}


template<class ParticleType>
const Foam::labelListList&
Foam::InteractionLists<ParticleType>::rilInverse() const
{
    return rilInverse_;
}


template<class ParticleType>
const Foam::labelListList& Foam::InteractionLists<ParticleType>::rwfil() const
{
    return rwfil_;
}


template<class ParticleType>
const Foam::labelListList&
Foam::InteractionLists<ParticleType>::rwfilInverse() const
{
    return rwfilInverse_;
}


template<class ParticleType>
const Foam::List<Foam::labelPair>&
Foam::InteractionLists<ParticleType>::cellIndexAndTransformToDistribute() const
{
    return cellIndexAndTransformToDistribute_;
}


template<class ParticleType>
const Foam::List<Foam::labelPair>&
Foam::InteractionLists<ParticleType>::
wallFaceIndexAndTransformToDistribute() const
{
    return wallFaceIndexAndTransformToDistribute_;
}


template<class ParticleType>
const Foam::List<Foam::referredWallFace>&
Foam::InteractionLists<ParticleType>::referredWallFaces() const
{
    return referredWallFaces_;
}


template<class ParticleType>
const Foam::word& Foam::InteractionLists<ParticleType>::UName() const
{
    return UName_;
}


template<class ParticleType>
const Foam::List<Foam::vector>&
Foam::InteractionLists<ParticleType>::referredWallData() const
{
    return referredWallData_;
}


template<class ParticleType>
const Foam::List<Foam::IDLList<ParticleType>>&
Foam::InteractionLists<ParticleType>::referredParticles() const
{
    return referredParticles_;
}


template<class ParticleType>
Foam::List<Foam::IDLList<ParticleType>>&
Foam::InteractionLists<ParticleType>::referredParticles()
{
    return referredParticles_;
}


// ************************************************************************* //
