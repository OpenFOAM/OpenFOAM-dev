/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2026 OpenFOAM Foundation
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

#include "NSRDS14Thermo.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class EquationOfState>
inline Foam::NSRDS14Thermo<EquationOfState>::NSRDS14Thermo
(
    const word& name,
    const NSRDS14Thermo& ns
)
:
    EquationOfState(name, ns),
    hf_(ns.hf_),
    sf_(ns.sf_),
    Tc_(ns.Tc_),
    a_(ns.a_),
    b_(ns.b_),
    c_(ns.c_),
    d_(ns.d_),
    p_({a_*a_, b_, -2*a_*c_, -a_*d_, -c_*c_/3, -0.5*c_*d_, -0.2*d_*d_}),
    hsRef_(ns.hsRef_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class EquationOfState>
inline Foam::scalar Foam::NSRDS14Thermo<EquationOfState>::limit
(
    const scalar T
) const
{
    return T;
}


template<class EquationOfState>
inline Foam::scalar Foam::NSRDS14Thermo<EquationOfState>::Cp
(
    const scalar p,
    const scalar T
) const
{
    return p_.value(max(1 - T/Tc_, small)) + EquationOfState::Cp(p, T);
}


template<class EquationOfState>
inline Foam::scalar Foam::NSRDS14Thermo<EquationOfState>::hs
(
    const scalar p,
    const scalar T
) const
{
    return
      - (p_.integral(max(1 - T/Tc_, small)) - hsRef_)*Tc_
      + EquationOfState::h(p, T);
}


template<class EquationOfState>
inline Foam::scalar Foam::NSRDS14Thermo<EquationOfState>::ha
(
    const scalar p,
    const scalar T
) const
{
    return hs(p, T) + hf();
}


template<class EquationOfState>
inline Foam::scalar Foam::NSRDS14Thermo<EquationOfState>::hf()
const
{
    return hf_;
}


template<class EquationOfState>
inline Foam::scalar Foam::NSRDS14Thermo<EquationOfState>::s
(
    const scalar p,
    const scalar T
) const
{
    NotImplemented;
    return NaN;
}


template<class EquationOfState>
inline Foam::scalar Foam::NSRDS14Thermo<EquationOfState>::gStd
(
    const scalar T
) const
{
    NotImplemented;
    return NaN;
}


template<class EquationOfState>
inline Foam::scalar Foam::NSRDS14Thermo<EquationOfState>::dCpdT
(
    const scalar p,
    const scalar T
) const
{
    return - p_.derivative(max(1 - T/Tc_, small))/Tc_;
}


// ************************************************************************* //
