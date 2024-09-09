/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2024 OpenFOAM Foundation
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

#include "chemistryModel.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ThermoType>
inline const Foam::multicomponentMixture<ThermoType>&
Foam::chemistryModel<ThermoType>::mixture() const
{
    return mixture_;
}


template<class ThermoType>
inline const Foam::PtrList<Foam::Reaction<ThermoType>>&
Foam::chemistryModel<ThermoType>::reactions() const
{
    return reactions_;
}


template<class ThermoType>
inline const Foam::PtrList<ThermoType>&
Foam::chemistryModel<ThermoType>::specieThermos() const
{
    return specieThermos_;
}


template<class ThermoType>
inline Foam::label Foam::chemistryModel<ThermoType>::nReaction() const
{
    return reactions_.size();
}


template<class ThermoType>
inline const Foam::word& Foam::chemistryModel<ThermoType>::reactionName
(
    const label reactioni
) const
{
    return reactions_[reactioni].name();
}


template<class ThermoType>
inline const Foam::PtrList<Foam::DimensionedField<Foam::scalar, Foam::volMesh>>&
Foam::chemistryModel<ThermoType>::RR() const
{
    return RR_;
}


// ************************************************************************* //
