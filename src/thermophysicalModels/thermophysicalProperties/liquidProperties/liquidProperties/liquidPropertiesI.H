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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

inline Foam::scalar Foam::liquidProperties::limit(const scalar T) const
{
    return T;
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


inline Foam::scalar Foam::liquidProperties::es(scalar p, const scalar T) const
{
    return hs(p, T);
}


inline Foam::scalar Foam::liquidProperties::ea(scalar p, const scalar T) const
{
    return ha(p, T);
}


template<class Func>
inline void Foam::liquidProperties::readIfPresent
(
    Func& f,
    const word& name,
    const dictionary& dict
)
{
    if (dict.found(name))
    {
        f = Func(name, dict.subDict(name));
    }
}


template<class Liquid>
inline void Foam::liquidProperties::readIfPresent
(
    Liquid& l,
    const dictionary& dict
)
{
    l.liquidProperties::readIfPresent(dict);
    readIfPresent(l.rho_, "rho", dict);
    readIfPresent(l.pv_, "pv", dict);
    readIfPresent(l.hl_, "hl", dict);
    readIfPresent(l.Cp_, "Cp", dict);
    readIfPresent(l.h_, "h", dict);
    readIfPresent(l.Cpg_, "Cpg", dict);
    readIfPresent(l.B_, "B", dict);
    readIfPresent(l.mu_, "mu", dict);
    readIfPresent(l.mug_, "mug", dict);
    readIfPresent(l.kappa_, "kappa", dict);
    readIfPresent(l.kappag_, "kappag", dict);
    readIfPresent(l.sigma_, "sigma", dict);
    // readIfPresent(l.D_, "D", dict);
}


template<class Liquid>
inline void Foam::liquidProperties::write
(
    const Liquid& l,
    Ostream& os
) const
{
    l.liquidProperties::write(os);
    writeEntry(os, l.rho_);
    writeEntry(os, l.pv_);
    writeEntry(os, l.hl_);
    writeEntry(os, l.Cp_);
    writeEntry(os, l.h_);
    writeEntry(os, l.Cpg_);
    writeEntry(os, l.B_);
    writeEntry(os, l.mu_);
    writeEntry(os, l.mug_);
    writeEntry(os,  l.kappa_);
    writeEntry(os, l.kappag_);
    writeEntry(os, l.sigma_);
    writeEntry(os,  l.D_);
}


// ************************************************************************* //
