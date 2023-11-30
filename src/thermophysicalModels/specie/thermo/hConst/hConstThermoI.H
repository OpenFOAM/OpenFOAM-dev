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

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class EquationOfState>
inline Foam::hConstThermo<EquationOfState>::hConstThermo
(
    const EquationOfState& st,
    const scalar Cp,
    const scalar hf,
    const scalar Tref,
    const scalar Href
)
:
    EquationOfState(st),
    Cp_(Cp),
    hf_(hf),
    Tref_(Tref),
    hsRef_(Href)
{}


template<class EquationOfState>
inline Foam::hConstThermo<EquationOfState>::hConstThermo
(
    const word& name,
    const hConstThermo& ct
)
:
    EquationOfState(name, ct),
    Cp_(ct.Cp_),
    hf_(ct.hf_),
    Tref_(ct.Tref_),
    hsRef_(ct.hsRef_)
{}


template<class EquationOfState>
inline Foam::autoPtr<Foam::hConstThermo<EquationOfState>>
Foam::hConstThermo<EquationOfState>::clone() const
{
    return autoPtr<hConstThermo<EquationOfState>>
    (
        new hConstThermo<EquationOfState>(*this)
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class EquationOfState>
inline Foam::scalar Foam::hConstThermo<EquationOfState>::limit
(
    const scalar T
) const
{
    return T;
}


template<class EquationOfState>
inline Foam::scalar Foam::hConstThermo<EquationOfState>::Cp
(
    const scalar p,
    const scalar T
) const
{
    return Cp_ + EquationOfState::Cp(p, T);
}


template<class EquationOfState>
inline Foam::scalar Foam::hConstThermo<EquationOfState>::hs
(
    const scalar p,
    const scalar T
) const
{
    return Cp_*(T - Tref_) + hsRef_ + EquationOfState::h(p, T);
}


template<class EquationOfState>
inline Foam::scalar Foam::hConstThermo<EquationOfState>::ha
(
    const scalar p,
    const scalar T
) const
{
    return hs(p, T) + hf();
}


template<class EquationOfState>
inline Foam::scalar Foam::hConstThermo<EquationOfState>::hf() const
{
    return hf_;
}


template<class EquationOfState>
inline Foam::scalar Foam::hConstThermo<EquationOfState>::s
(
    const scalar p,
    const scalar T
) const
{
    return Cp_*log(T/Tstd) + EquationOfState::sp(p, T);
}


template<class EquationOfState>
inline Foam::scalar Foam::hConstThermo<EquationOfState>::gStd
(
    const scalar T
) const
{
    return Cp_*(T - Tref_) + hsRef_ + hf() - Cp_*T*log(T/Tstd);
}


template<class EquationOfState>
inline Foam::scalar Foam::hConstThermo<EquationOfState>::dCpdT
(
    const scalar p,
    const scalar T
) const
{
    return 0;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class EquationOfState>
inline void Foam::hConstThermo<EquationOfState>::operator+=
(
    const hConstThermo<EquationOfState>& ct
)
{
    scalar Y1 = this->Y();

    EquationOfState::operator+=(ct);

    if (mag(this->Y()) > small)
    {
        if
        (
            hConstThermo<EquationOfState>::debug
         && notEqual(Tref_, ct.Tref_)
        )
        {
            FatalErrorInFunction
                << "Tref " << Tref_ << " for "
                << (this->name().size() ? this->name() : "others")
                << " != " << ct.Tref_ << " for "
                << (ct.name().size() ? ct.name() : "others")
                << exit(FatalError);
        }

        Y1 /= this->Y();
        const scalar Y2 = ct.Y()/this->Y();

        Cp_ = Y1*Cp_ + Y2*ct.Cp_;
        hf_ = Y1*hf_ + Y2*ct.hf_;
        hsRef_ = Y1*hsRef_ + Y2*ct.hsRef_;
    }
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template<class EquationOfState>
inline Foam::hConstThermo<EquationOfState> Foam::operator+
(
    const hConstThermo<EquationOfState>& ct1,
    const hConstThermo<EquationOfState>& ct2
)
{
    EquationOfState eofs
    (
        static_cast<const EquationOfState&>(ct1)
      + static_cast<const EquationOfState&>(ct2)
    );

    if (mag(eofs.Y()) < small)
    {
        return hConstThermo<EquationOfState>
        (
            eofs,
            ct1.Cp_,
            ct1.hf_,
            ct1.Tref_,
            ct1.hsRef_
        );
    }
    else
    {
        if
        (
            hConstThermo<EquationOfState>::debug
         && notEqual(ct1.Tref_, ct2.Tref_)
        )
        {
            FatalErrorInFunction
                << "Tref " << ct1.Tref_ << " for "
                << (ct1.name().size() ? ct1.name() : "others")
                << " != " << ct2.Tref_ << " for "
                << (ct2.name().size() ? ct2.name() : "others")
                << exit(FatalError);
        }

        return hConstThermo<EquationOfState>
        (
            eofs,
            ct1.Y()/eofs.Y()*ct1.Cp_
          + ct2.Y()/eofs.Y()*ct2.Cp_,
            ct1.Y()/eofs.Y()*ct1.hf_
          + ct2.Y()/eofs.Y()*ct2.hf_,
            ct1.Tref_,
            ct1.Y()/eofs.Y()*ct1.hsRef_
          + ct2.Y()/eofs.Y()*ct2.hsRef_
        );
    }
}


template<class EquationOfState>
inline Foam::hConstThermo<EquationOfState> Foam::operator*
(
    const scalar s,
    const hConstThermo<EquationOfState>& ct
)
{
    return hConstThermo<EquationOfState>
    (
        s*static_cast<const EquationOfState&>(ct),
        ct.Cp_,
        ct.hf_,
        ct.Tref_,
        ct.hsRef_
    );
}


template<class EquationOfState>
inline Foam::hConstThermo<EquationOfState> Foam::operator==
(
    const hConstThermo<EquationOfState>& ct1,
    const hConstThermo<EquationOfState>& ct2
)
{
    EquationOfState eofs
    (
        static_cast<const EquationOfState&>(ct1)
     == static_cast<const EquationOfState&>(ct2)
    );

    if
    (
        hConstThermo<EquationOfState>::debug
     && notEqual(ct1.Tref_, ct2.Tref_)
    )
    {
        FatalErrorInFunction
            << "Tref " << ct1.Tref_ << " for "
            << (ct1.name().size() ? ct1.name() : "others")
            << " != " << ct2.Tref_ << " for "
            << (ct2.name().size() ? ct2.name() : "others")
            << exit(FatalError);
    }

    return hConstThermo<EquationOfState>
    (
        eofs,
        ct2.Y()/eofs.Y()*ct2.Cp_
      - ct1.Y()/eofs.Y()*ct1.Cp_,
        ct2.Y()/eofs.Y()*ct2.hf_
      - ct1.Y()/eofs.Y()*ct1.hf_,
        ct1.Tref_,
        ct2.Y()/eofs.Y()*ct2.hsRef_
      - ct1.Y()/eofs.Y()*ct1.hsRef_
    );
}


// ************************************************************************* //
