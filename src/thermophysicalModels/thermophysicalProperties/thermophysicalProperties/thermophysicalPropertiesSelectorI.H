/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2017-2023 OpenFOAM Foundation
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

#include "thermophysicalPropertiesSelector.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ThermophysicalProperties>
inline const Foam::word&
Foam::thermophysicalPropertiesSelector<ThermophysicalProperties>::name() const
{
    return name_;
}


template<class ThermophysicalProperties>
inline const ThermophysicalProperties&
Foam::thermophysicalPropertiesSelector<ThermophysicalProperties>::
properties() const
{
    return propertiesPtr_();
}


template<class ThermophysicalProperties>
inline Foam::scalar
Foam::thermophysicalPropertiesSelector<ThermophysicalProperties>::W() const
{
    return propertiesPtr_->W();
}


template<class ThermophysicalProperties>
inline Foam::scalar
Foam::thermophysicalPropertiesSelector<ThermophysicalProperties>::limit
(
    const scalar T
) const
{
    return propertiesPtr_->limit(T);
}


template<class ThermophysicalProperties>
inline Foam::scalar
Foam::thermophysicalPropertiesSelector<ThermophysicalProperties>::rho
(
    scalar p,
    scalar T
) const
{
    return propertiesPtr_->rho(p, T);
}


template<class ThermophysicalProperties>
inline Foam::scalar
Foam::thermophysicalPropertiesSelector<ThermophysicalProperties>::psi
(
    scalar p,
    scalar T
) const
{
    return propertiesPtr_->psi(p, T);
}


template<class ThermophysicalProperties>
inline Foam::scalar
Foam::thermophysicalPropertiesSelector<ThermophysicalProperties>::CpMCv
(
    scalar p,
    scalar T
) const
{
    return propertiesPtr_->CpMCv(p, T);
}


template<class ThermophysicalProperties>
inline Foam::scalar
Foam::thermophysicalPropertiesSelector<ThermophysicalProperties>::Cp
(
    scalar p,
    scalar T
) const
{
    return propertiesPtr_->Cp(p, T);
}


template<class ThermophysicalProperties>
inline Foam::scalar
Foam::thermophysicalPropertiesSelector<ThermophysicalProperties>::Hs
(
    scalar p,
    scalar T
) const
{
    return propertiesPtr_->Hs(p, T);
}


template<class ThermophysicalProperties>
inline Foam::scalar
Foam::thermophysicalPropertiesSelector<ThermophysicalProperties>::Hf() const
{
    return propertiesPtr_->Hf();
}


template<class ThermophysicalProperties>
inline Foam::scalar
Foam::thermophysicalPropertiesSelector<ThermophysicalProperties>::Ha
(
    scalar p,
    scalar T
) const
{
    return propertiesPtr_->Ha(p, T);
}


template<class ThermophysicalProperties>
inline Foam::scalar
Foam::thermophysicalPropertiesSelector<ThermophysicalProperties>::Cv
(
    scalar p,
    scalar T
) const
{
    return propertiesPtr_->Cp(p, T) - CpMCv(p, T);
}


template<class ThermophysicalProperties>
inline Foam::scalar
Foam::thermophysicalPropertiesSelector<ThermophysicalProperties>::Es
(
    scalar p,
    scalar T
) const
{
    return propertiesPtr_->Es(p, T);
}


template<class ThermophysicalProperties>
inline Foam::scalar
Foam::thermophysicalPropertiesSelector<ThermophysicalProperties>::Ea
(
    scalar p,
    scalar T
) const
{
    return propertiesPtr_->Ea(p, T);
}


template<class ThermophysicalProperties>
inline Foam::scalar
Foam::thermophysicalPropertiesSelector<ThermophysicalProperties>::S
(
    scalar p,
    scalar T
) const
{
    return propertiesPtr_->S(p, T);
}


template<class ThermophysicalProperties>
inline Foam::scalar
Foam::thermophysicalPropertiesSelector<ThermophysicalProperties>::kappa
(
    scalar p,
    scalar T
) const
{
    return propertiesPtr_->kappa(p, T);
}


// ************************************************************************* //
