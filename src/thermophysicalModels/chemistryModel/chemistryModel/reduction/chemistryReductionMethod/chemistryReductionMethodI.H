/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2021 OpenFOAM Foundation
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

#include "chemistryReductionMethod.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ThermoType>
inline Foam::label Foam::chemistryReductionMethod<ThermoType>::nSpecie()
{
    return nSpecie_;
}


template<class ThermoType>
inline Foam::label
Foam::chemistryReductionMethod<ThermoType>::nActiveSpecies() const
{
    return nActiveSpecies_;
}


template<class ThermoType>
inline const Foam::List<bool>&
Foam::chemistryReductionMethod<ThermoType>::activeSpecies() const
{
    return activeSpecies_;
}


template<class ThermoType>
inline Foam::scalar
Foam::chemistryReductionMethod<ThermoType>::tolerance() const
{
    return tolerance_;
}


template<class ThermoType>
inline bool Foam::chemistryReductionMethod<ThermoType>::reactionDisabled
(
    const label i
) const
{
    return reactionsDisabled_[i];
}


// ************************************************************************* //
