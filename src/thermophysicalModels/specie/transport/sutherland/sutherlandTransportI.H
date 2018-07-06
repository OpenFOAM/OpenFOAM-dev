/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Thermo>
inline void Foam::sutherlandTransport<Thermo>::calcCoeffs
(
    const scalar mu1, const scalar T1,
    const scalar mu2, const scalar T2
)
{
    scalar rootT1 = sqrt(T1);
    scalar mu1rootT2 = mu1*sqrt(T2);
    scalar mu2rootT1 = mu2*rootT1;

    Ts_ = (mu2rootT1 - mu1rootT2)/(mu1rootT2/T1 - mu2rootT1/T2);

    As_ = mu1*(1.0 + Ts_/T1)/rootT1;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Thermo>
inline Foam::sutherlandTransport<Thermo>::sutherlandTransport
(
    const Thermo& t,
    const scalar As,
    const scalar Ts
)
:
    Thermo(t),
    As_(As),
    Ts_(Ts)
{}


template<class Thermo>
inline Foam::sutherlandTransport<Thermo>::sutherlandTransport
(
    const Thermo& t,
    const scalar mu1, const scalar T1,
    const scalar mu2, const scalar T2
)
:
    Thermo(t)
{
    calcCoeffs(mu1, T1, mu2, T2);
}


template<class Thermo>
inline Foam::sutherlandTransport<Thermo>::sutherlandTransport
(
    const word& name,
    const sutherlandTransport& st
)
:
    Thermo(name, st),
    As_(st.As_),
    Ts_(st.Ts_)
{}


template<class Thermo>
inline Foam::autoPtr<Foam::sutherlandTransport<Thermo>>
Foam::sutherlandTransport<Thermo>::clone() const
{
    return autoPtr<sutherlandTransport<Thermo>>
    (
        new sutherlandTransport<Thermo>(*this)
    );
}


template<class Thermo>
inline Foam::autoPtr<Foam::sutherlandTransport<Thermo>>
Foam::sutherlandTransport<Thermo>::New
(
    const dictionary& dict
)
{
    return autoPtr<sutherlandTransport<Thermo>>
    (
        new sutherlandTransport<Thermo>(dict)
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Thermo>
inline Foam::scalar Foam::sutherlandTransport<Thermo>::mu
(
    const scalar p,
    const scalar T
) const
{
    return As_*::sqrt(T)/(1.0 + Ts_/T);
}


template<class Thermo>
inline Foam::scalar Foam::sutherlandTransport<Thermo>::kappa
(
    const scalar p, const scalar T
) const
{
    scalar Cv_ = this->Cv(p, T);
    return mu(p, T)*Cv_*(1.32 + 1.77*this->R()/Cv_);
}


template<class Thermo>
inline Foam::scalar Foam::sutherlandTransport<Thermo>::alphah
(
    const scalar p,
    const scalar T
) const
{

    return kappa(p, T)/this->Cp(p, T);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Thermo>
inline void Foam::sutherlandTransport<Thermo>::operator=
(
    const sutherlandTransport<Thermo>& st
)
{
    Thermo::operator=(st);

    As_ = st.As_;
    Ts_ = st.Ts_;
}


template<class Thermo>
inline void Foam::sutherlandTransport<Thermo>::operator+=
(
    const sutherlandTransport<Thermo>& st
)
{
    scalar Y1 = this->Y();

    Thermo::operator+=(st);

    if (mag(this->Y()) > small)
    {
        Y1 /= this->Y();
        scalar Y2 = st.Y()/this->Y();

        As_ = Y1*As_ + Y2*st.As_;
        Ts_ = Y1*Ts_ + Y2*st.Ts_;
    }
}


template<class Thermo>
inline void Foam::sutherlandTransport<Thermo>::operator*=
(
    const scalar s
)
{
    Thermo::operator*=(s);
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template<class Thermo>
inline Foam::sutherlandTransport<Thermo> Foam::operator+
(
    const sutherlandTransport<Thermo>& st1,
    const sutherlandTransport<Thermo>& st2
)
{
    Thermo t
    (
        static_cast<const Thermo&>(st1) + static_cast<const Thermo&>(st2)
    );

    if (mag(t.Y()) < small)
    {
        return sutherlandTransport<Thermo>
        (
            t,
            0,
            st1.As_,
            st1.Ts_
        );
    }
    else
    {
        scalar Y1 = st1.Y()/t.Y();
        scalar Y2 = st2.Y()/t.Y();

        return sutherlandTransport<Thermo>
        (
            t,
            Y1*st1.As_ + Y2*st2.As_,
            Y1*st1.Ts_ + Y2*st2.Ts_
        );
    }
}


template<class Thermo>
inline Foam::sutherlandTransport<Thermo> Foam::operator*
(
    const scalar s,
    const sutherlandTransport<Thermo>& st
)
{
    return sutherlandTransport<Thermo>
    (
        s*static_cast<const Thermo&>(st),
        st.As_,
        st.Ts_
    );
}


// ************************************************************************* //
