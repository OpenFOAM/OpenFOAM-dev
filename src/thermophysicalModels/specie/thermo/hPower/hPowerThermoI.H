/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2012-2023 OpenFOAM Foundation
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

#include "hPowerThermo.H"
#include "specie.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class EquationOfState>
inline void Foam::hPowerThermo<EquationOfState>::checkT
(
    const scalar T
) const
{
    if (T < 0)
    {
        FatalErrorInFunction
            << "attempt to evaluate hPowerThermo<EquationOfState>"
               " for negative temperature " << T
            << abort(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class EquationOfState>
inline Foam::hPowerThermo<EquationOfState>::hPowerThermo
(
    const EquationOfState& st,
    const scalar c0,
    const scalar n0,
    const scalar Tref,
    const scalar hf
)
:
    EquationOfState(st),
    c0_(c0),
    n0_(n0),
    Tref_(Tref),
    hf_(hf)
{}


template<class EquationOfState>
inline Foam::hPowerThermo<EquationOfState>::hPowerThermo
(
    const word& name,
    const hPowerThermo& jt
)
:
    EquationOfState(name, jt),
    c0_(jt.c0_),
    n0_(jt.n0_),
    Tref_(jt.Tref_),
    hf_(jt.hf_)
{}


template<class EquationOfState>
inline Foam::autoPtr<Foam::hPowerThermo<EquationOfState>>
Foam::hPowerThermo<EquationOfState>::clone() const
{
    return autoPtr<hPowerThermo<EquationOfState>>
    (
        new hPowerThermo<EquationOfState>(*this)
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class EquationOfState>
inline Foam::scalar Foam::hPowerThermo<EquationOfState>::limit
(
    const scalar T
) const
{
    return T;
}


template<class EquationOfState>
inline Foam::scalar Foam::hPowerThermo<EquationOfState>::Cp
(
    const scalar p,
    const scalar T
) const
{
    return c0_*pow(T/Tref_, n0_) + EquationOfState::Cp(p, T);
}


template<class EquationOfState>
inline Foam::scalar Foam::hPowerThermo<EquationOfState>::hs
(
    const scalar p,
    const scalar T
) const
{
    return
        c0_*(pow(T, n0_ + 1) - pow(Tstd, n0_ + 1))/(pow(Tref_, n0_)*(n0_ + 1))
      + EquationOfState::h(p, T);
}


template<class EquationOfState>
inline Foam::scalar Foam::hPowerThermo<EquationOfState>::ha
(
    const scalar p,
    const scalar T
) const
{
    return hs(p, T) + hf();
}


template<class EquationOfState>
inline Foam::scalar Foam::hPowerThermo<EquationOfState>::hf() const
{
    return hf_;
}


template<class EquationOfState>
inline Foam::scalar Foam::hPowerThermo<EquationOfState>::s
(
    const scalar p,
    const scalar T
) const
{
    return
        c0_*(pow(T, n0_) - pow(Tstd, n0_))/(pow(Tref_, n0_)*n0_)
      + EquationOfState::sp(p, T);
}


template<class EquationOfState>
inline Foam::scalar Foam::hPowerThermo<EquationOfState>::gStd
(
    const scalar T
) const
{
    return
        c0_*(pow(T, n0_ + 1) - pow(Tstd, n0_ + 1))/(pow(Tref_, n0_)*(n0_ + 1))
      + hf()
      - c0_*(pow(T, n0_) - pow(Tstd, n0_))*T/(pow(Tref_, n0_)*n0_);
}


template<class EquationOfState>
inline Foam::scalar Foam::hPowerThermo<EquationOfState>::dCpdT
(
    const scalar p,
    const scalar T
) const
{
    // To be implemented
    NotImplemented;
    return 0;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class EquationOfState>
inline void Foam::hPowerThermo<EquationOfState>::operator+=
(
    const hPowerThermo<EquationOfState>& ct
)
{
    scalar Y1 = this->Y();

    EquationOfState::operator+=(ct);

    if (mag(this->Y()) > small)
    {
        Y1 /= this->Y();
        const scalar Y2 = ct.Y()/this->Y();

        hf_ = Y1*hf_ + Y2*ct.hf_;
        c0_ = Y1*c0_ + Y2*ct.c0_;
        n0_ = Y1*n0_ + Y2*ct.n0_;
        Tref_ = Y1*Tref_ + Y2*ct.Tref_;
    }
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template<class EquationOfState>
inline Foam::hPowerThermo<EquationOfState> Foam::operator+
(
    const hPowerThermo<EquationOfState>& ct1,
    const hPowerThermo<EquationOfState>& ct2
)
{
    EquationOfState eofs
    (
        static_cast<const EquationOfState&>(ct1)
      + static_cast<const EquationOfState&>(ct2)
    );

    if (mag(eofs.Y()) < small)
    {
        return hPowerThermo<EquationOfState>
        (
            eofs,
            ct1.c0_,
            ct1.n0_,
            ct1.Tref_,
            ct1.hf_
        );
    }
    else
    {
        return hPowerThermo<EquationOfState>
        (
            eofs,
            ct1.Y()/eofs.Y()*ct1.c0_
          + ct2.Y()/eofs.Y()*ct2.c0_,
            ct1.Y()/eofs.Y()*ct1.n0_
          + ct2.Y()/eofs.Y()*ct2.n0_,
            ct1.Y()/eofs.Y()*ct1.Tref_
          + ct2.Y()/eofs.Y()*ct2.Tref_,
            ct1.Y()/eofs.Y()*ct1.hf_
          + ct2.Y()/eofs.Y()*ct2.hf_
        );
    }
}


template<class EquationOfState>
inline Foam::hPowerThermo<EquationOfState> Foam::operator*
(
    const scalar s,
    const hPowerThermo<EquationOfState>& ct
)
{
    return hPowerThermo<EquationOfState>
    (
        s*static_cast<const EquationOfState&>(ct),
        ct.c0_,
        ct.n0_,
        ct.Tref_,
        ct.hf_
    );
}


template<class EquationOfState>
inline Foam::hPowerThermo<EquationOfState> Foam::operator==
(
    const hPowerThermo<EquationOfState>& ct1,
    const hPowerThermo<EquationOfState>& ct2
)
{
    EquationOfState eofs
    (
        static_cast<const EquationOfState&>(ct1)
     == static_cast<const EquationOfState&>(ct2)
    );

    return hPowerThermo<EquationOfState>
    (
        eofs,
        ct2.Y()/eofs.Y()*ct2.c0_
      - ct1.Y()/eofs.Y()*ct1.c0_,
        ct2.Y()/eofs.Y()*ct2.n0_
      - ct1.Y()/eofs.Y()*ct1.n0_,
        ct2.Y()/eofs.Y()*ct2.Tref_
      - ct1.Y()/eofs.Y()*ct1.Tref_,
        ct2.Y()/eofs.Y()*ct2.hf_
      - ct1.Y()/eofs.Y()*ct1.hf_
    );
}


// ************************************************************************* //
