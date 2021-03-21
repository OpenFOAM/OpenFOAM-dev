/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2021 OpenFOAM Foundation
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

inline Foam::scalar Foam::C14H30::rho(scalar p, scalar T) const
{
    return rho_.value(T);
}


inline Foam::scalar Foam::C14H30::pv(scalar p, scalar T) const
{
    return pv_.value(T);
}


inline Foam::scalar Foam::C14H30::hl(scalar p, scalar T) const
{
    return hl_.value(T);
}


inline Foam::scalar Foam::C14H30::Cp(scalar p, scalar T) const
{
    return Cp_.value(T);
}


inline Foam::scalar Foam::C14H30::Hs(scalar p, scalar T) const
{
    return Ha(p, T) - Hf();
}


inline Foam::scalar Foam::C14H30::Hf() const
{
    return Hf_;
}


inline Foam::scalar Foam::C14H30::Ha(scalar p, scalar T) const
{
    return h_.value(T);
}


inline Foam::scalar Foam::C14H30::Cpg(scalar p, scalar T) const
{
    return Cpg_.value(T);
}


inline Foam::scalar Foam::C14H30::B(scalar p, scalar T) const
{
    return B_.value(T);
}


inline Foam::scalar Foam::C14H30::mu(scalar p, scalar T) const
{
    return mu_.value(T);
}


inline Foam::scalar Foam::C14H30::mug(scalar p, scalar T) const
{
    return mug_.value(T);
}


inline Foam::scalar Foam::C14H30::kappa(scalar p, scalar T) const
{
    return kappa_.value(T);
}


inline Foam::scalar Foam::C14H30::kappag(scalar p, scalar T) const
{
    return kappag_.value(T);
}


inline Foam::scalar Foam::C14H30::sigma(scalar p, scalar T) const
{
    return sigma_.value(T);
}


inline Foam::scalar Foam::C14H30::D(scalar p, scalar T) const
{
    return D_.value(p, T);
}


inline Foam::scalar Foam::C14H30::D(scalar p, scalar T, scalar Wb) const
{
    return D_.value(p, T, Wb);
}


// ************************************************************************* //
