/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2017 OpenFOAM Foundation
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

inline Foam::scalar Foam::liquidProperties::limit(const scalar T) const
{
    return T;
}


inline Foam::scalar Foam::liquidProperties::W() const
{
    return W_;
}


inline Foam::scalar Foam::liquidProperties::Y() const
{
    return 1;
}


inline Foam::scalar Foam::liquidProperties::Tc() const
{
    return Tc_;
}


inline Foam::scalar Foam::liquidProperties::Pc() const
{
    return Pc_;
}


inline Foam::scalar Foam::liquidProperties::Vc() const
{
    return Vc_;
}


inline Foam::scalar Foam::liquidProperties::Zc() const
{
    return Zc_;
}


inline Foam::scalar Foam::liquidProperties::Tt() const
{
    return Tt_;
}


inline Foam::scalar Foam::liquidProperties::Pt() const
{
    return Pt_;
}


inline Foam::scalar Foam::liquidProperties::Tb() const
{
    return Tb_;
}


inline Foam::scalar Foam::liquidProperties::dipm() const
{
    return dipm_;
}


inline Foam::scalar Foam::liquidProperties::omega() const
{
    return omega_;
}


inline Foam::scalar Foam::liquidProperties::delta() const
{
    return delta_;
}


inline Foam::scalar Foam::liquidProperties::psi(scalar p, scalar T) const
{
    return 0;
}


inline Foam::scalar Foam::liquidProperties::CpMCv(scalar p, scalar T) const
{
    return 0;
}


inline Foam::scalar Foam::liquidProperties::Ha(scalar p, scalar T) const
{
    return h(p, T);
}


inline Foam::scalar Foam::liquidProperties::Hs(scalar p, scalar T) const
{
    return h(p, T);
}


inline Foam::scalar Foam::liquidProperties::Hc() const
{
    return 0;
}


inline Foam::scalar Foam::liquidProperties::alphah(scalar p, scalar T) const
{
    return kappa(p, T)/Cp(p, T);
}


// ************************************************************************* //
