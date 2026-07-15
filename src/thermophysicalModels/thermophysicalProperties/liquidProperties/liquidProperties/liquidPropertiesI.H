/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2026 OpenFOAM Foundation
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

#include "liquidProperties.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

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


inline Foam::scalar Foam::liquidProperties::limit(const scalar T) const
{
    return T;
}


inline Foam::scalar
Foam::liquidProperties::psi(const scalar p, const scalar T) const
{
    return 0;
}


inline Foam::scalar
Foam::liquidProperties::CpMCv(const scalar p, const scalar T) const
{
    return 0;
}


inline Foam::scalar
Foam::liquidProperties::es(const scalar p, const scalar T) const
{
    return hs(p, T);
}


inline Foam::scalar
Foam::liquidProperties::ea(const scalar p, const scalar T) const
{
    return ha(p, T);
}


template<class Func>
inline void Foam::liquidProperties::readIfPresent
(
    Func& f,
    const Function1<scalar>&,
    const word& name,
    const dimensionSet& dims,
    const dictionary& dict
)
{
    if (dict.found(name))
    {
        f = Func(name, {dimensions::temperature, dims}, dict.subDict(name));
    }
}


template<class Func>
inline void Foam::liquidProperties::readIfPresent
(
    Func& f,
    const Function2<scalar>&,
    const word& name,
    const dimensionSet& dims,
    const dictionary& dict
)
{
    if (dict.found(name))
    {
        f = Func
        (
            name,
            {dimensions::temperature, dimensions::pressure, dims},
            dict.subDict(name)
        );
    }
}


template<class Func>
inline void Foam::liquidProperties::readIfPresent
(
    Func& f,
    const word& name,
    const dimensionSet& dims,
    const dictionary& dict
)
{
    readIfPresent(f, f, name, dims, dict);
}


template<class Liquid>
inline void Foam::liquidProperties::readIfPresent
(
    Liquid& l,
    const dictionary& dict
)
{
    l.liquidProperties::readIfPresent(dict);
    readIfPresent(l.rho_, "rho", dimensions::density, dict);
    readIfPresent(l.pv_, "pv", dimensions::pressure, dict);
    readIfPresent(l.hl_, "hl", dimensions::specificEnergy, dict);
    readIfPresent(l.Cp_, "Cp", dimensions::specificHeatCapacity, dict);
    readIfPresent(l.h_, "h", dimensions::specificEnergy, dict);
    readIfPresent(l.Cpg_, "Cpg", dimensions::specificHeatCapacity, dict);
    readIfPresent(l.B_, "B", dimensions::volume/dimensions::mass, dict);
    readIfPresent(l.mu_, "mu", dimensions::dynamicViscosity, dict);
    readIfPresent(l.mug_, "mug", dimensions::dynamicViscosity, dict);
    readIfPresent(l.kappa_, "kappa", dimensions::thermalConductivity, dict);
    readIfPresent(l.kappag_, "kappag", dimensions::thermalConductivity, dict);
    readIfPresent
    (
        l.sigma_,
        "sigma",
        dimensions::force/dimensions::length,
        dict
    );
    readIfPresent(l.D_, "D", dimensions::kinematicDiffusivity, dict);
}


template<class Func>
inline void Foam::liquidProperties::write
(
    const Func& f,
    const Function1<scalar>&,
    const dimensionSet& dims,
    Ostream& os
) const
{
    writeEntry(os, dimensions::temperature, dims, f);
}


template<class Func>
inline void Foam::liquidProperties::write
(
    const Func& f,
    const Function2<scalar>&,
    const dimensionSet& dims,
    Ostream& os
) const
{
    writeEntry(os, dimensions::pressure, dimensions::temperature, dims, f);
}


template<class Func>
inline void Foam::liquidProperties::write
(
    const Func& f,
    const dimensionSet& dims,
    Ostream& os
) const
{
    write(f, f, dims, os);
}


template<class Func>
inline void Foam::liquidProperties::write
(
    const autoPtr<Func>& fPtr,
    const dimensionSet& dims,
    Ostream& os
) const
{
    write(fPtr(), fPtr(), dims, os);
}


template<class Liquid>
inline void Foam::liquidProperties::write
(
    const Liquid& l,
    Ostream& os
) const
{
    l.liquidProperties::write(os);
    write(l.rho_, dimensions::density, os);
    write(l.pv_, dimensions::pressure, os);
    write(l.hl_, dimensions::specificEnergy, os);
    write(l.Cp_, dimensions::specificHeatCapacity, os);
    write(l.h_, dimensions::specificEnergy, os);
    write(l.Cpg_, dimensions::specificHeatCapacity, os);
    write(l.B_, dimensions::volume/dimensions::mass, os);
    write(l.mu_, dimensions::dynamicViscosity, os);
    write(l.mug_, dimensions::dynamicViscosity, os);
    write(l.kappa_, dimensions::thermalConductivity, os);
    write(l.kappag_, dimensions::thermalConductivity, os);
    write(l.sigma_, dimensions::force/dimensions::length, os);
    write(l.D_, dimensions::kinematicDiffusivity, os);
}


// ************************************************************************* //
