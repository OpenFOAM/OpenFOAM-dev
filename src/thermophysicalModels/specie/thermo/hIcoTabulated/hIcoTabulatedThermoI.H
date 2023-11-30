/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2020-2023 OpenFOAM Foundation
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

#include "hIcoTabulatedThermo.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class EquationOfState>
inline Foam::hIcoTabulatedThermo<EquationOfState>::hIcoTabulatedThermo
(
    const word& name,
    const hIcoTabulatedThermo& pt
)
:
    EquationOfState(name, pt),
    hf_(pt.hf_),
    sf_(pt.sf_),
    Cp_(pt.Cp_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class EquationOfState>
inline Foam::scalar Foam::hIcoTabulatedThermo<EquationOfState>::limit
(
    const scalar T
) const
{
    return T;
}


template<class EquationOfState>
inline Foam::scalar Foam::hIcoTabulatedThermo<EquationOfState>::Cp
(
    const scalar p,
    const scalar T
) const
{
    return Cp_.value(T) + EquationOfState::Cp(p, T);
}


template<class EquationOfState>
inline Foam::scalar Foam::hIcoTabulatedThermo<EquationOfState>::hs
(
    const scalar p,
    const scalar T
) const
{
    return Cp_.intfdT(T) + EquationOfState::h(p, T);
}


template<class EquationOfState>
inline Foam::scalar Foam::hIcoTabulatedThermo<EquationOfState>::ha
(
    const scalar p,
    const scalar T
) const
{
    return hs(p, T) + hf_;
}


template<class EquationOfState>
inline Foam::scalar Foam::hIcoTabulatedThermo<EquationOfState>::hf()
const
{
    return hf_;
}


template<class EquationOfState>
inline Foam::scalar Foam::hIcoTabulatedThermo<EquationOfState>::s
(
    const scalar p,
    const scalar T
) const
{
    return Cp_.intfByTdT(T) + EquationOfState::sp(p, T) + sf_;
}


template<class EquationOfState>
inline Foam::scalar Foam::hIcoTabulatedThermo<EquationOfState>::gStd
(
    const scalar T
) const
{
    return Cp_.intfdT(T) + hf_ - (Cp_.intfByTdT(T) + sf_)*T;
}


template<class EquationOfState>
inline Foam::scalar Foam::hIcoTabulatedThermo<EquationOfState>::dCpdT
(
    const scalar p,
    const scalar T
) const
{
    return Cp_.dfdT(T);
}


// ************************************************************************* //
