/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2015-2018 OpenFOAM Foundation
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

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class EquationOfState>
inline Foam::hRefConstThermo<EquationOfState>::hRefConstThermo
(
    const EquationOfState& st,
    const scalar cp,
    const scalar hf,
    const scalar tref,
    const scalar href
)
:
    EquationOfState(st),
    Cp_(cp),
    Hf_(hf),
    Tref_(tref),
    Href_(href)
{}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class EquationOfState>
inline Foam::hRefConstThermo<EquationOfState>::hRefConstThermo
(
    const word& name,
    const hRefConstThermo& ct
)
:
    EquationOfState(name, ct),
    Cp_(ct.Cp_),
    Hf_(ct.Hf_),
    Tref_(ct.Tref_),
    Href_(ct.Href_)
{}


template<class EquationOfState>
inline Foam::autoPtr<Foam::hRefConstThermo<EquationOfState>>
Foam::hRefConstThermo<EquationOfState>::clone() const
{
    return autoPtr<hRefConstThermo<EquationOfState>>
    (
        new hRefConstThermo<EquationOfState>(*this)
    );
}


template<class EquationOfState>
inline Foam::autoPtr<Foam::hRefConstThermo<EquationOfState>>
Foam::hRefConstThermo<EquationOfState>::New(const dictionary& dict)
{
    return autoPtr<hRefConstThermo<EquationOfState>>
    (
        new hRefConstThermo<EquationOfState>(dict)
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class EquationOfState>
inline Foam::scalar Foam::hRefConstThermo<EquationOfState>::limit
(
    const scalar T
) const
{
    return T;
}


template<class EquationOfState>
inline Foam::scalar Foam::hRefConstThermo<EquationOfState>::Cp
(
    const scalar p,
    const scalar T
) const
{
    return Cp_ + EquationOfState::Cp(p, T);
}


template<class EquationOfState>
inline Foam::scalar Foam::hRefConstThermo<EquationOfState>::Ha
(
    const scalar p, const scalar T
) const
{
    return Cp_*(T-Tref_) + Href_ + Hf_ + EquationOfState::H(p, T);
}


template<class EquationOfState>
inline Foam::scalar Foam::hRefConstThermo<EquationOfState>::Hs
(
    const scalar p, const scalar T
) const
{
    return Cp_*(T-Tref_) + Href_ + EquationOfState::H(p, T);
}


template<class EquationOfState>
inline Foam::scalar Foam::hRefConstThermo<EquationOfState>::Hc() const
{
    return Hf_;
}


template<class EquationOfState>
inline Foam::scalar Foam::hRefConstThermo<EquationOfState>::S
(
    const scalar p, const scalar T
) const
{
    return Cp_*log(T/Tstd) + EquationOfState::S(p, T);
}


template<class EquationOfState>
inline Foam::scalar Foam::hRefConstThermo<EquationOfState>::dGdT
(
    const scalar p, const scalar T
) const
{
    return 0;
}


template<class EquationOfState>
inline Foam::scalar Foam::hRefConstThermo<EquationOfState>::dCpdT
(
    const scalar p, const scalar T
) const
{
    return 0;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class EquationOfState>
inline void Foam::hRefConstThermo<EquationOfState>::operator+=
(
    const hRefConstThermo<EquationOfState>& ct
)
{
    scalar Y1 = this->Y();

    EquationOfState::operator+=(ct);

    if (mag(this->Y()) > small)
    {
        Y1 /= this->Y();
        const scalar Y2 = ct.Y()/this->Y();

        Cp_ = Y1*Cp_ + Y2*ct.Cp_;
        Hf_ = Y1*Hf_ + Y2*ct.Hf_;
    }
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template<class EquationOfState>
inline Foam::hRefConstThermo<EquationOfState> Foam::operator+
(
    const hRefConstThermo<EquationOfState>& ct1,
    const hRefConstThermo<EquationOfState>& ct2
)
{
    EquationOfState eofs
    (
        static_cast<const EquationOfState&>(ct1)
      + static_cast<const EquationOfState&>(ct2)
    );

    if (mag(eofs.Y()) < small)
    {
        return hRefConstThermo<EquationOfState>
        (
            eofs,
            ct1.Cp_,
            ct1.Hf_,
            ct1.Tref_,
            ct1.Href_
        );
    }
    else
    {
        return hRefConstThermo<EquationOfState>
        (
            eofs,
            ct1.Y()/eofs.Y()*ct1.Cp_
          + ct2.Y()/eofs.Y()*ct2.Cp_,
            ct1.Y()/eofs.Y()*ct1.Hf_
          + ct2.Y()/eofs.Y()*ct2.Hf_,
            ct1.Y()/eofs.Y()*ct1.Tref_
          + ct2.Y()/eofs.Y()*ct2.Tref_,
            ct1.Y()/eofs.Y()*ct1.Href_
          + ct2.Y()/eofs.Y()*ct2.Href_
        );
    }
}


template<class EquationOfState>
inline Foam::hRefConstThermo<EquationOfState> Foam::operator*
(
    const scalar s,
    const hRefConstThermo<EquationOfState>& ct
)
{
    return hRefConstThermo<EquationOfState>
    (
        s*static_cast<const EquationOfState&>(ct),
        ct.Cp_,
        ct.Hf_,
        ct.Tref_,
        ct.Href_
    );
}


template<class EquationOfState>
inline Foam::hRefConstThermo<EquationOfState> Foam::operator==
(
    const hRefConstThermo<EquationOfState>& ct1,
    const hRefConstThermo<EquationOfState>& ct2
)
{
    EquationOfState eofs
    (
        static_cast<const EquationOfState&>(ct1)
     == static_cast<const EquationOfState&>(ct2)
    );

    return hRefConstThermo<EquationOfState>
    (
        eofs,
        ct2.Y()/eofs.Y()*ct2.Cp_
      - ct1.Y()/eofs.Y()*ct1.Cp_,
        ct2.Y()/eofs.Y()*ct2.Hf_
      - ct1.Y()/eofs.Y()*ct1.Hf_
    );
}


// ************************************************************************* //
