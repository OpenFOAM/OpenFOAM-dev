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

#include "eIcoTabulatedThermo.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class EquationOfState>
inline Foam::eIcoTabulatedThermo<EquationOfState>::eIcoTabulatedThermo
(
    const word& name,
    const eIcoTabulatedThermo& pt
)
:
    EquationOfState(name, pt),
    hf_(pt.hf_),
    sf_(pt.sf_),
    Cv_(pt.Cv_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class EquationOfState>
inline Foam::scalar Foam::eIcoTabulatedThermo<EquationOfState>::limit
(
    const scalar T
) const
{
    return T;
}


template<class EquationOfState>
inline Foam::scalar Foam::eIcoTabulatedThermo<EquationOfState>::Cv
(
    const scalar p,
    const scalar T
) const
{
    return Cv_.value(T) + EquationOfState::Cv(p, T);
}


template<class EquationOfState>
inline Foam::scalar Foam::eIcoTabulatedThermo<EquationOfState>::es
(
    const scalar p,
    const scalar T
) const
{
    return Cv_.intfdT(T) + EquationOfState::e(p, T);
}


template<class EquationOfState>
inline Foam::scalar Foam::eIcoTabulatedThermo<EquationOfState>::ea
(
    const scalar p,
    const scalar T
) const
{
    return es(p, T) + hf_;
}


template<class EquationOfState>
inline Foam::scalar Foam::eIcoTabulatedThermo<EquationOfState>::hf()
const
{
    return hf_;
}


template<class EquationOfState>
inline Foam::scalar Foam::eIcoTabulatedThermo<EquationOfState>::s
(
    const scalar p,
    const scalar T
) const
{
    NotImplemented;
    return 0;
}


template<class EquationOfState>
inline Foam::scalar Foam::eIcoTabulatedThermo<EquationOfState>::gStd
(
    const scalar T
) const
{
    NotImplemented;
    return 0;
}


template<class EquationOfState>
inline Foam::scalar Foam::eIcoTabulatedThermo<EquationOfState>::dCpdT
(
    const scalar p,
    const scalar T
) const
{
    NotImplemented;
    return 0;
}


// ************************************************************************* //
