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

#include "NSRDSThermo.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class EquationOfState>
inline Foam::NSRDSThermo<EquationOfState>::NSRDSThermo
(
    const EquationOfState& ns,
    const scalar hf,
    const scalar sf,
    const FixedLaurentPolynomial<scalar, 0, 5>& CpCoeffs,
    const scalar hsRef,
    const scalar sRef
)
:
    EquationOfState(ns),
    hf_(hf),
    sf_(sf),
    CpCoeffs_(CpCoeffs),
    hsRef_(hsRef),
    sRef_(sRef)
{}


template<class EquationOfState>
inline Foam::NSRDSThermo<EquationOfState>::NSRDSThermo
(
    const word& name,
    const NSRDSThermo& ns
)
:
    EquationOfState(name, ns),
    hf_(ns.hf_),
    sf_(ns.sf_),
    CpCoeffs_(ns.CpCoeffs_),
    hsRef_(ns.hsRef_),
    sRef_(ns.hsRef_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class EquationOfState>
inline Foam::scalar Foam::NSRDSThermo<EquationOfState>::limit
(
    const scalar T
) const
{
    return T;
}


template<class EquationOfState>
inline Foam::scalar Foam::NSRDSThermo<EquationOfState>::Cp
(
    const scalar p,
    const scalar T
) const
{
    return CpCoeffs_.value(T) + EquationOfState::Cp(p, T);
}


template<class EquationOfState>
inline Foam::scalar Foam::NSRDSThermo<EquationOfState>::hs
(
    const scalar p,
    const scalar T
) const
{
    return CpCoeffs_.integral(T) - hsRef_ + EquationOfState::h(p, T);
}


template<class EquationOfState>
inline Foam::scalar Foam::NSRDSThermo<EquationOfState>::ha
(
    const scalar p,
    const scalar T
) const
{
    return hs(p, T) + hf();
}


template<class EquationOfState>
inline Foam::scalar Foam::NSRDSThermo<EquationOfState>::hf()
const
{
    return hf_;
}


template<class EquationOfState>
inline Foam::scalar Foam::NSRDSThermo<EquationOfState>::s
(
    const scalar p,
    const scalar T
) const
{
    return CpCoeffs_.byX().integral(T) - sRef_ + EquationOfState::sp(p, T);
}


template<class EquationOfState>
inline Foam::scalar Foam::NSRDSThermo<EquationOfState>::gStd
(
    const scalar T
) const
{
    return
        CpCoeffs_.integral(T) - hsRef_ + hf()
      - (CpCoeffs_.byX().integral(T) - sRef_)*T;
}


template<class EquationOfState>
inline Foam::scalar Foam::NSRDSThermo<EquationOfState>::dCpdT
(
    const scalar p,
    const scalar T
) const
{
    return CpCoeffs_.derivative(T);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class EquationOfState>
inline void Foam::NSRDSThermo<EquationOfState>::operator+=
(
    const NSRDSThermo<EquationOfState>& ns
)
{
    scalar Y1 = this->Y();

    EquationOfState::operator+=(ns);

    if (mag(this->Y()) > small)
    {
        Y1 /= this->Y();
        const scalar Y2 = ns.Y()/this->Y();

        hf_ = Y1*hf_ + Y2*ns.hf_;
        sf_ = Y1*sf_ + Y2*ns.sf_;
        CpCoeffs_ = Y1*CpCoeffs_ + Y2*ns.CpCoeffs_;
        hsRef_ = Y1*hsRef_ + Y2*ns.hsRef_;
        sRef_ = Y1*sRef_ + Y2*ns.sRef_;
    }
}


template<class EquationOfState>
inline void Foam::NSRDSThermo<EquationOfState>::operator*=
(
    const scalar s
)
{
    EquationOfState::operator*=(s);
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template<class EquationOfState>
inline Foam::NSRDSThermo<EquationOfState> Foam::operator+
(
    const NSRDSThermo<EquationOfState>& ns1,
    const NSRDSThermo<EquationOfState>& ns2
)
{
    EquationOfState eofs = ns1;
    eofs += ns2;

    if (mag(eofs.Y()) < small)
    {
        return NSRDSThermo<EquationOfState>
        (
            eofs,
            ns1.hf_,
            ns1.sf_,
            ns1.CpCoeffs_,
            ns1.hsRef_,
            ns1.sRef_
        );
    }
    {
        const scalar Y1 = ns1.Y()/eofs.Y();
        const scalar Y2 = ns2.Y()/eofs.Y();

        return NSRDSThermo<EquationOfState>
        (
            eofs,
            Y1*ns1.hf_ + Y2*ns2.hf_,
            Y1*ns1.sf_ + Y2*ns2.sf_,
            Y1*ns1.CpCoeffs_ + Y2*ns2.CpCoeffs_,
            Y1*ns1.hsRef_ + Y2*ns2.hsRef_,
            Y1*ns1.sRef_ + Y2*ns2.sRef_
        );
    }
}


template<class EquationOfState>
inline Foam::NSRDSThermo<EquationOfState> Foam::operator*
(
    const scalar s,
    const NSRDSThermo<EquationOfState>& ns
)
{
    return NSRDSThermo<EquationOfState>
    (
        s*static_cast<const EquationOfState&>(ns),
        ns.hf_,
        ns.sf_,
        ns.CpCoeffs_,
        ns.hsRef_,
        ns.sRef_
    );
}


template<class EquationOfState>
inline Foam::NSRDSThermo<EquationOfState> Foam::operator==
(
    const NSRDSThermo<EquationOfState>& ns1,
    const NSRDSThermo<EquationOfState>& ns2
)
{
    EquationOfState eofs
    (
        static_cast<const EquationOfState&>(ns1)
     == static_cast<const EquationOfState&>(ns2)
    );

    const scalar Y1 = ns1.Y()/eofs.Y();
    const scalar Y2 = ns2.Y()/eofs.Y();

    return NSRDSThermo<EquationOfState>
    (
        eofs,
        Y2*ns2.hf_       - Y1*ns1.hf_,
        Y2*ns2.sf_       - Y1*ns1.sf_,
        Y2*ns2.CpCoeffs_ - Y1*ns1.CpCoeffs_,
        Y2*ns2.hsRef_    - Y1*ns1.hsRef_,
        Y2*ns2.sRef_     - Y1*ns1.sRef_
    );
}


// ************************************************************************* //
