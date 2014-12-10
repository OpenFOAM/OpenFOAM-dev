/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2014  OpenFOAM Foundation
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

#include "PengRobinsonGas.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Specie>
inline Foam::PengRobinsonGas<Specie>::PengRobinsonGas
(
    const Specie& sp,
    const scalar& Tc,
    const scalar& Vc,
    const scalar& Zc,
    const scalar& Pc,
    const scalar& omega
)
:
    Specie(sp),
    Tc_(Tc),
    Vc_(Vc),
    Zc_(Zc),
    Pc_(Pc),
    omega_(omega)
{}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Specie>
inline Foam::PengRobinsonGas<Specie>::PengRobinsonGas
(
    const word& name,
    const PengRobinsonGas& pg
)
:
    Specie(name, pg),
    Tc_(pg.Tc_),
    Pc_(pg.Pc_),
    Vc_(pg.Vc_),
    Zc_(pg.Zc_),
    omega_(pg.omega_)
{}


template<class Specie>
inline Foam::autoPtr<Foam::PengRobinsonGas <Specie> >
Foam::PengRobinsonGas<Specie>::clone() const
{
    return autoPtr<PengRobinsonGas<Specie> >
    (
        new PengRobinsonGas<Specie>(*this)
    );
}


template<class Specie>
inline Foam::autoPtr<Foam::PengRobinsonGas<Specie> >
Foam::PengRobinsonGas<Specie>::New
(
    Istream& is
)
{
    return autoPtr<PengRobinsonGas<Specie> >(new PengRobinsonGas<Specie>(is));
}


template<class Specie>
inline Foam::autoPtr<Foam::PengRobinsonGas<Specie> >
Foam::PengRobinsonGas<Specie>::New
(
    const dictionary& dict
)
{
    return autoPtr<PengRobinsonGas<Specie> >
    (
        new PengRobinsonGas<Specie>(dict)
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Specie>
inline Foam::scalar Foam::PengRobinsonGas<Specie>::rho
(
    scalar p,
    scalar T
) const
{
    scalar z = Z(p, T);
    return p/(z*this->R()*T);
}


template<class Specie>
inline Foam::scalar Foam::PengRobinsonGas<Specie>::psi
(
    scalar p,
    scalar T
) const
{
    scalar z = Z(p, T);
    return 1.0/(z*this->R()*T);
}


template<class Specie>
inline Foam::scalar Foam::PengRobinsonGas<Specie>::Z
(
    scalar p,
    scalar T
) const
{
    scalar a = 0.45724*sqr(this->R())*sqr(Tc_)/Pc_;
    scalar b = 0.07780*this->R()*Tc_/Pc_;
    scalar Tr = T/Tc_;

    scalar alpha =
        sqr
        (
            1.0
          + (0.37464 + 1.54226*omega_- 0.26992*sqr(omega_))
          * (1.0 - sqrt(Tr))
        );

    scalar B = b*p/(this->R()*T);
    scalar A = a*alpha*p/sqr(this->R()*T);

    scalar a2 = B - 1.0;
    scalar a1 = A - 2.0*B - 3.0*sqr(B);
    scalar a0 = -A*B + sqr(B) + pow3(B);

    scalar Q = (3.0*a1 - a2*a2)/9.0;
    scalar Rl = (9.0*a2*a1 - 27.0*a0 - 2.0*a2*a2*a2)/54;

    scalar Q3 = Q*Q*Q;
    scalar D = Q3 + Rl*Rl;

    scalar root = -1.0;

    if (D <= 0.0)
    {
        scalar th = ::acos(Rl/sqrt(-Q3));
        scalar qm = 2.0*sqrt(-Q);
        scalar r1 = qm*cos(th/3.0) - a2/3.0;
        scalar r2 = qm*cos((th + 2.0*constant::mathematical::pi)/3.0) - a2/3.0;
        scalar r3 = qm*cos((th + 4.0*constant::mathematical::pi)/3.0) - a2/3.0;

        root = max(r1, max(r2, r3));
    }
    else
    {
        // one root is real
        scalar D05 = sqrt(D);
        scalar S = pow(Rl + D05, 1.0/3.0);
        scalar Tl = 0;
        if (D05 > Rl)
        {
            Tl = -pow(mag(Rl - D05), 1.0/3.0);
        }
        else
        {
            Tl = pow(Rl - D05, 1.0/3.0);
        }

        root = S + Tl - a2/3.0;
    }

    return root;
}


template<class Specie>
inline Foam::scalar Foam::PengRobinsonGas<Specie>::cpMcv
(
    scalar p,
    scalar T
) const
{
    return this->RR*Z(p, T);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Specie>
inline void Foam::PengRobinsonGas<Specie>::operator+=
(
    const PengRobinsonGas<Specie>& pg
)
{
    scalar molr1 = this->nMoles();
    Specie::operator+=(pg);

    molr1 /= this->nMoles();
    scalar molr2 = pg.nMoles()/this->nMoles();

    Tc_ = molr1*Tc_ + molr2*pg.Tc_;
    Vc_ = molr1*Vc_ + molr2*pg.Vc_;
    Zc_ = molr1*Zc_ + molr2*pg.Zc_;
    Pc_ = specie::RR*Zc_*Tc_/Vc_;
    omega_ = molr1*omega_ + molr2*pg.omega_;
}


template<class Specie>
inline void Foam::PengRobinsonGas<Specie>::operator-=
(
    const PengRobinsonGas<Specie>& pg
)
{
    scalar molr1 = this->nMoles();

    Specie::operator-=(pg);

    molr1 /= this->nMoles();
    scalar molr2 = pg.nMoles()/this->nMoles();

    Tc_ = molr1*Tc_ - molr2*pg.Tc_;
    Vc_ = molr1*Vc_ - molr2*pg.Vc_;
    Zc_ = molr1*Zc_ - molr2*pg.Zc_;
    Pc_ = specie::RR*Zc_*Tc_/Vc_;
    omega_ = molr1*omega_ - molr2*pg.omega_;
}


template<class Specie>
inline void Foam::PengRobinsonGas<Specie>::operator*=(const scalar s)
{
     Specie::operator*=(s);
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //


template<class Specie>
Foam::PengRobinsonGas<Specie> Foam::operator+
(
    const PengRobinsonGas<Specie>& pg1,
    const PengRobinsonGas<Specie>& pg2
)
{
    scalar nMoles = pg1.nMoles() + pg2.nMoles();
    scalar molr1 = pg1.nMoles()/nMoles;
    scalar molr2 = pg2.nMoles()/nMoles;

    const scalar Tc = molr1*pg1.Tc_ + molr2*pg2.Tc_;
    const scalar Vc = molr1*pg1.Vc_ + molr2*pg2.Vc_;
    const scalar Zc = molr1*pg1.Zc_ + molr2*pg2.Zc_;

    return PengRobinsonGas<Specie>
    (
        static_cast<const Specie&>(pg1)
      + static_cast<const Specie&>(pg2),
        Tc,
        Vc,
        Zc,
        specie::RR*Zc*Tc/Vc,
        molr1*pg1.omega_ + molr2*pg2.omega_
    );
}


template<class Specie>
Foam::PengRobinsonGas<Specie> Foam::operator-
(
    const PengRobinsonGas<Specie>& pg1,
    const PengRobinsonGas<Specie>& pg2
)
{
    scalar nMoles = pg1.nMoles() + pg2.nMoles();
    scalar molr1 = pg1.nMoles()/nMoles;
    scalar molr2 = pg2.nMoles()/nMoles;

    const scalar Tc = molr1*pg1.Tc_ + molr2*pg2.Tc_;
    const scalar Vc = molr1*pg1.Vc_ + molr2*pg2.Vc_;
    const scalar Zc = molr1*pg1.Zc_ + molr2*pg2.Zc_;

    return PengRobinsonGas<Specie>
    (
        static_cast<const Specie&>(pg1)
      - static_cast<const Specie&>(pg2),
        Tc,
        Vc,
        Zc,
        specie::RR*Zc*Tc/Vc,
        molr1*pg1.omega_ - molr2*pg2.omega_
    );
}


template<class Specie>
Foam::PengRobinsonGas<Specie> Foam::operator*
(
    const scalar s,
    const PengRobinsonGas<Specie>& pg
)
{
    return PengRobinsonGas<Specie>
    (
        s*static_cast<const Specie&>(pg),
        pg.Tc_,
        pg.Vc_,
        pg.Zc_,
        pg.Pc_,
        pg.omega_
    );
}


template<class Specie>
Foam::PengRobinsonGas<Specie> Foam::operator==
(
    const PengRobinsonGas<Specie>& pg1,
    const PengRobinsonGas<Specie>& pg2
)
{
    return pg2 - pg1;
}


// ************************************************************************* //
