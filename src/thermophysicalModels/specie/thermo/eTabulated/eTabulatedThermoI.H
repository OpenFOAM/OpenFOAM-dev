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

#include "eTabulatedThermo.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class EquationOfState>
inline Foam::eTabulatedThermo<EquationOfState>::eTabulatedThermo
(
    const word& name,
    const eTabulatedThermo& pt
)
:
    EquationOfState(name, pt),
    hf_(pt.hf_),
    sf_(pt.sf_),
    es_(pt.es_),
    Cp_(pt.Cp_),
    Cv_(pt.Cv_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class EquationOfState>
inline Foam::scalar Foam::eTabulatedThermo<EquationOfState>::limit
(
    const scalar T
) const
{
    return T;
}


template<class EquationOfState>
inline Foam::scalar Foam::eTabulatedThermo<EquationOfState>::Cp
(
    const scalar p,
    const scalar T
) const
{
    return Cp_.value(p, T);
}


template<class EquationOfState>
inline Foam::scalar Foam::eTabulatedThermo<EquationOfState>::Cv
(
    const scalar p,
    const scalar T
) const
{
    return Cv_.value(p, T);
}


template<class EquationOfState>
inline Foam::scalar Foam::eTabulatedThermo<EquationOfState>::es
(
    const scalar p,
    const scalar T
) const
{
    return es_.value(p, T);
}


template<class EquationOfState>
inline Foam::scalar Foam::eTabulatedThermo<EquationOfState>::ea
(
    const scalar p,
    const scalar T
) const
{
    return es(p, T) + hf_;
}


template<class EquationOfState>
inline Foam::scalar Foam::eTabulatedThermo<EquationOfState>::hs
(
    const scalar p,
    const scalar T
) const
{
    return es(p, T) + p/EquationOfState::rho(p, T);
}


template<class EquationOfState>
inline Foam::scalar Foam::eTabulatedThermo<EquationOfState>::ha
(
    const scalar p,
    const scalar T
) const
{
    return ea(p, T) + p/EquationOfState::rho(p, T);
}


template<class EquationOfState>
inline Foam::scalar Foam::eTabulatedThermo<EquationOfState>::hf()
const
{
    return hf_;
}


template<class EquationOfState>
inline Foam::scalar Foam::eTabulatedThermo<EquationOfState>::s
(
    const scalar p,
    const scalar T
) const
{
    NotImplemented;
    return 0;
}


template<class EquationOfState>
inline Foam::scalar Foam::eTabulatedThermo<EquationOfState>::gStd
(
    const scalar T
) const
{
    NotImplemented;
    return 0;
}


template<class EquationOfState>
inline Foam::scalar Foam::eTabulatedThermo<EquationOfState>::dCpdT
(
    const scalar p,
    const scalar T
) const
{
    NotImplemented;
    return 0;
}


// ************************************************************************* //
