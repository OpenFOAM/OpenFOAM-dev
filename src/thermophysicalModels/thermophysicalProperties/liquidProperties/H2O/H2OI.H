/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2023 OpenFOAM Foundation
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

inline Foam::scalar Foam::H2O::rho(scalar p, scalar T) const
{
    return rho_.value(T);
}


inline Foam::scalar Foam::H2O::pv(scalar p, scalar T) const
{
    return pv_.value(T);
}


inline Foam::scalar Foam::H2O::hl(scalar p, scalar T) const
{
    return hl_.value(T);
}


inline Foam::scalar Foam::H2O::Cp(scalar p, scalar T) const
{
    return Cp_.value(T);
}


inline Foam::scalar Foam::H2O::hs(scalar p, scalar T) const
{
    return ha(p, T) - hf();
}


inline Foam::scalar Foam::H2O::hf() const
{
    return hf_;
}


inline Foam::scalar Foam::H2O::ha(scalar p, scalar T) const
{
    return h_.value(T);
}


inline Foam::scalar Foam::H2O::Cpg(scalar p, scalar T) const
{
    return Cpg_.value(T);
}


inline Foam::scalar Foam::H2O::B(scalar p, scalar T) const
{
    return B_.value(T);
}


inline Foam::scalar Foam::H2O::mu(scalar p, scalar T) const
{
    return mu_.value(T);
}


inline Foam::scalar Foam::H2O::mug(scalar p, scalar T) const
{
    return mug_.value(T);
}


inline Foam::scalar Foam::H2O::kappa(scalar p, scalar T) const
{
    return kappa_.value(T);
}


inline Foam::scalar Foam::H2O::kappag(scalar p, scalar T) const
{
    return kappag_.value(T);
}


inline Foam::scalar Foam::H2O::sigma(scalar p, scalar T) const
{
    return sigma_.value(T);
}


inline Foam::scalar Foam::H2O::D(scalar p, scalar T) const
{
    return D_.value(p, T);
}


inline Foam::scalar Foam::H2O::D(scalar p, scalar T, scalar Wb) const
{
    return D_.value(p, T, Wb);
}


// ************************************************************************* //
