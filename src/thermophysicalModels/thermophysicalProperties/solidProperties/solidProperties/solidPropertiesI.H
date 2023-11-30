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

#include "specie.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

inline Foam::scalar Foam::solidProperties::rho() const
{
    return rho_;
}


inline Foam::scalar Foam::solidProperties::Cp() const
{
    return Cp_;
}


inline Foam::scalar Foam::solidProperties::hf() const
{
    return hf_;
}


inline Foam::scalar Foam::solidProperties::hs(const scalar T) const
{
    return Cp_*(T - Tstd);
}


inline Foam::scalar Foam::solidProperties::ha(const scalar T) const
{
    return hs(T) + hf();
}


inline Foam::scalar Foam::solidProperties::es(const scalar T) const
{
    return hs(T);
}


inline Foam::scalar Foam::solidProperties::ea(const scalar T) const
{
    return es(T) + hf();
}


inline Foam::scalar Foam::solidProperties::kappa() const
{
    return kappa_;
}


inline Foam::scalar Foam::solidProperties::emissivity() const
{
    return emissivity_;
}


// ************************************************************************* //
