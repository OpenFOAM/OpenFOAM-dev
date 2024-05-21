/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2024 OpenFOAM Foundation
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
    const Function1<scalar>&,
    const word& name,
    const dimensionSet& dims,
    const dictionary& dict
)
{
    if (dict.found(name))
    {
        f = Func(name, {dimTemperature, dims}, dict.subDict(name));
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
        f = Func(name, {dimTemperature, dimPressure, dims}, dict.subDict(name));
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
    readIfPresent(l.rho_, "rho", dimDensity, dict);
    readIfPresent(l.pv_, "pv", dimPressure, dict);
    readIfPresent(l.hl_, "hl", dimEnergy/dimMass, dict);
    readIfPresent(l.Cp_, "Cp", dimSpecificHeatCapacity, dict);
    readIfPresent(l.h_, "h", dimEnergy/dimMass, dict);
    readIfPresent(l.Cpg_, "Cpg", dimSpecificHeatCapacity, dict);
    readIfPresent(l.B_, "B", dimVolume/dimMass, dict);
    readIfPresent(l.mu_, "mu", dimDynamicViscosity, dict);
    readIfPresent(l.mug_, "mug", dimDynamicViscosity, dict);
    readIfPresent(l.kappa_, "kappa", dimThermalConductivity, dict);
    readIfPresent(l.kappag_, "kappag", dimThermalConductivity, dict);
    readIfPresent(l.sigma_, "sigma", dimForce/dimLength, dict);
    readIfPresent(l.D_, "D", dimArea/dimTime, dict);
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
    writeEntry(os, dimTemperature, dims, f);
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
    writeEntry(os, dimPressure, dimTemperature, dims, f);
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
    write(l.rho_, dimDensity, os);
    write(l.pv_, dimPressure, os);
    write(l.hl_, dimEnergy/dimMass, os);
    write(l.Cp_, dimSpecificHeatCapacity, os);
    write(l.h_, dimEnergy/dimMass, os);
    write(l.Cpg_, dimSpecificHeatCapacity, os);
    write(l.B_, dimVolume/dimMass, os);
    write(l.mu_, dimDynamicViscosity, os);
    write(l.mug_, dimDynamicViscosity, os);
    write(l.kappa_, dimThermalConductivity, os);
    write(l.kappag_, dimThermalConductivity, os);
    write(l.sigma_, dimForce/dimLength, os);
    write(l.D_, dimArea/dimTime, os);
}


// ************************************************************************* //
