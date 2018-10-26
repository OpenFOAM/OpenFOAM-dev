/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2018 OpenFOAM Foundation
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
inline Foam::eRefConstThermo<EquationOfState>::eRefConstThermo
(
    const EquationOfState& st,
    const scalar cv,
    const scalar hf,
    const scalar tref,
    const scalar eref
)
:
    EquationOfState(st),
    Cv_(cv),
    Hf_(hf),
    Tref_(tref),
    Eref_(eref)
{}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class EquationOfState>
inline Foam::eRefConstThermo<EquationOfState>::eRefConstThermo
(
    const word& name,
    const eRefConstThermo& ct
)
:
    EquationOfState(name, ct),
    Cv_(ct.Cv_),
    Hf_(ct.Hf_),
    Tref_(ct.Tref_),
    Eref_(ct.Eref_)
{}


template<class EquationOfState>
inline Foam::autoPtr<Foam::eRefConstThermo<EquationOfState>>
Foam::eRefConstThermo<EquationOfState>::clone() const
{
    return autoPtr<eRefConstThermo<EquationOfState>>
    (
        new eRefConstThermo<EquationOfState>(*this)
    );
}


template<class EquationOfState>
inline Foam::autoPtr<Foam::eRefConstThermo<EquationOfState>>
Foam::eRefConstThermo<EquationOfState>::New(const dictionary& dict)
{
    return autoPtr<eRefConstThermo<EquationOfState>>
    (
        new eRefConstThermo<EquationOfState>(dict)
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class EquationOfState>
inline Foam::scalar Foam::eRefConstThermo<EquationOfState>::limit
(
    const scalar T
) const
{
    return T;
}


template<class EquationOfState>
inline Foam::scalar Foam::eRefConstThermo<EquationOfState>::Cv
(
    const scalar p,
    const scalar T
) const
{
    return Cv_ + EquationOfState::Cv(p, T);
}


template<class EquationOfState>
inline Foam::scalar Foam::eRefConstThermo<EquationOfState>::Es
(
    const scalar p, const scalar T
) const
{
    return Cv_*(T - Tref_) + Eref_ + EquationOfState::E(p, T);
}


template<class EquationOfState>
inline Foam::scalar Foam::eRefConstThermo<EquationOfState>::Hc() const
{
    return Hf_;
}


template<class EquationOfState>
inline Foam::scalar Foam::eRefConstThermo<EquationOfState>::Ea
(
    const scalar p, const scalar T
) const
{
    return Es(p, T) + Hc();
}


template<class EquationOfState>
inline Foam::scalar Foam::eRefConstThermo<EquationOfState>::S
(
    const scalar p, const scalar T
) const
{
    return Cp(p, T)*log(T/Tstd) + EquationOfState::S(p, T);
}


template<class EquationOfState>
inline Foam::scalar Foam::eRefConstThermo<EquationOfState>::dGdT
(
    const scalar p, const scalar T
) const
{
    return 0;
}


template<class EquationOfState>
inline Foam::scalar Foam::eRefConstThermo<EquationOfState>::dCpdT
(
    const scalar p, const scalar T
) const
{
    return 0;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class EquationOfState>
inline void Foam::eRefConstThermo<EquationOfState>::operator+=
(
    const eRefConstThermo<EquationOfState>& ct
)
{
    scalar Y1 = this->Y();

    EquationOfState::operator+=(ct);

    if (mag(this->Y()) > small)
    {
        Y1 /= this->Y();
        const scalar Y2 = ct.Y()/this->Y();

        Cv_ = Y1*Cv_ + Y2*ct.Cv_;
        Hf_ = Y1*Hf_ + Y2*ct.Hf_;
    }
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template<class EquationOfState>
inline Foam::eRefConstThermo<EquationOfState> Foam::operator+
(
    const eRefConstThermo<EquationOfState>& ct1,
    const eRefConstThermo<EquationOfState>& ct2
)
{
    EquationOfState eofs
    (
        static_cast<const EquationOfState&>(ct1)
      + static_cast<const EquationOfState&>(ct2)
    );

    if (mag(eofs.Y()) < small)
    {
        return eRefConstThermo<EquationOfState>
        (
            eofs,
            ct1.Cv_,
            ct1.Hf_,
            ct1.Tref_,
            ct1.Eref_
        );
    }
    else
    {
        return eRefConstThermo<EquationOfState>
        (
            eofs,
            ct1.Y()/eofs.Y()*ct1.Cv_
          + ct2.Y()/eofs.Y()*ct2.Cv_,
            ct1.Y()/eofs.Y()*ct1.Hf_
          + ct2.Y()/eofs.Y()*ct2.Hf_,
            ct1.Y()/eofs.Y()*ct1.Tref_
          + ct2.Y()/eofs.Y()*ct2.Tref_,
            ct1.Y()/eofs.Y()*ct1.Eref_
          + ct2.Y()/eofs.Y()*ct2.Eref_
        );
    }
}


template<class EquationOfState>
inline Foam::eRefConstThermo<EquationOfState> Foam::operator*
(
    const scalar s,
    const eRefConstThermo<EquationOfState>& ct
)
{
    return eRefConstThermo<EquationOfState>
    (
        s*static_cast<const EquationOfState&>(ct),
        ct.Cv_,
        ct.Hf_,
        ct.Tref_,
        ct.Eref_
    );
}


template<class EquationOfState>
inline Foam::eRefConstThermo<EquationOfState> Foam::operator==
(
    const eRefConstThermo<EquationOfState>& ct1,
    const eRefConstThermo<EquationOfState>& ct2
)
{
    EquationOfState eofs
    (
        static_cast<const EquationOfState&>(ct1)
     == static_cast<const EquationOfState&>(ct2)
    );

    return eRefConstThermo<EquationOfState>
    (
        eofs,
        ct2.Y()/eofs.Y()*ct2.Cv_
      - ct1.Y()/eofs.Y()*ct1.Cv_,
        ct2.Y()/eofs.Y()*ct2.Hf_
      - ct1.Y()/eofs.Y()*ct1.Hf_
    );
}


// ************************************************************************* //
