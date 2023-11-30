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

#include "hPolynomialThermo.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class EquationOfState, int PolySize>
inline Foam::hPolynomialThermo<EquationOfState, PolySize>::hPolynomialThermo
(
    const EquationOfState& pt,
    const scalar hf,
    const scalar sf,
    const Polynomial<PolySize>& CpCoeffs,
    const typename Polynomial<PolySize>::intPolyType& hCoeffs,
    const Polynomial<PolySize>& sCoeffs
)
:
    EquationOfState(pt),
    hf_(hf),
    sf_(sf),
    CpCoeffs_(CpCoeffs),
    hCoeffs_(hCoeffs),
    sCoeffs_(sCoeffs)
{}


template<class EquationOfState, int PolySize>
inline Foam::hPolynomialThermo<EquationOfState, PolySize>::hPolynomialThermo
(
    const word& name,
    const hPolynomialThermo& pt
)
:
    EquationOfState(name, pt),
    hf_(pt.hf_),
    sf_(pt.sf_),
    CpCoeffs_(pt.CpCoeffs_),
    hCoeffs_(pt.hCoeffs_),
    sCoeffs_(pt.sCoeffs_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class EquationOfState, int PolySize>
inline Foam::scalar Foam::hPolynomialThermo<EquationOfState, PolySize>::limit
(
    const scalar T
) const
{
    return T;
}


template<class EquationOfState, int PolySize>
inline Foam::scalar Foam::hPolynomialThermo<EquationOfState, PolySize>::Cp
(
    const scalar p,
    const scalar T
) const
{
    return CpCoeffs_.value(T) + EquationOfState::Cp(p, T);
}


template<class EquationOfState, int PolySize>
inline Foam::scalar Foam::hPolynomialThermo<EquationOfState, PolySize>::hs
(
    const scalar p,
    const scalar T
) const
{
    return ha(p, T) - hf();
}


template<class EquationOfState, int PolySize>
inline Foam::scalar Foam::hPolynomialThermo<EquationOfState, PolySize>::ha
(
    const scalar p,
    const scalar T
) const
{
    return hCoeffs_.value(T) + EquationOfState::h(p, T);
}


template<class EquationOfState, int PolySize>
inline Foam::scalar Foam::hPolynomialThermo<EquationOfState, PolySize>::hf()
const
{
    return hf_;
}


template<class EquationOfState, int PolySize>
inline Foam::scalar Foam::hPolynomialThermo<EquationOfState, PolySize>::s
(
    const scalar p,
    const scalar T
) const
{
    return sCoeffs_.value(T) + EquationOfState::sp(p, T);
}


template<class EquationOfState, int PolySize>
inline Foam::scalar Foam::hPolynomialThermo<EquationOfState, PolySize>::gStd
(
    const scalar T
) const
{
    return hCoeffs_.value(T) - sCoeffs_.value(T)*T;
}


template<class EquationOfState, int PolySize>
inline Foam::scalar Foam::hPolynomialThermo<EquationOfState, PolySize>::dCpdT
(
    const scalar p,
    const scalar T
) const
{
    return
    (
        CpCoeffs_.derivative(T)
    );
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class EquationOfState, int PolySize>
inline void Foam::hPolynomialThermo<EquationOfState, PolySize>::operator+=
(
    const hPolynomialThermo<EquationOfState, PolySize>& pt
)
{
    scalar Y1 = this->Y();

    EquationOfState::operator+=(pt);

    if (mag(this->Y()) > small)
    {
        Y1 /= this->Y();
        const scalar Y2 = pt.Y()/this->Y();

        hf_ = Y1*hf_ + Y2*pt.hf_;
        sf_ = Y1*sf_ + Y2*pt.sf_;
        CpCoeffs_ = Y1*CpCoeffs_ + Y2*pt.CpCoeffs_;
        hCoeffs_ = Y1*hCoeffs_ + Y2*pt.hCoeffs_;
        sCoeffs_ = Y1*sCoeffs_ + Y2*pt.sCoeffs_;
    }
}


template<class EquationOfState, int PolySize>
inline void Foam::hPolynomialThermo<EquationOfState, PolySize>::operator*=
(
    const scalar s
)
{
    EquationOfState::operator*=(s);
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template<class EquationOfState, int PolySize>
inline Foam::hPolynomialThermo<EquationOfState, PolySize> Foam::operator+
(
    const hPolynomialThermo<EquationOfState, PolySize>& pt1,
    const hPolynomialThermo<EquationOfState, PolySize>& pt2
)
{
    EquationOfState eofs = pt1;
    eofs += pt2;

    if (mag(eofs.Y()) < small)
    {
        return hPolynomialThermo<EquationOfState, PolySize>
        (
            eofs,
            pt1.hf_,
            pt1.sf_,
            pt1.CpCoeffs_,
            pt1.hCoeffs_,
            pt1.sCoeffs_
        );
    }
    {
        const scalar Y1 = pt1.Y()/eofs.Y();
        const scalar Y2 = pt2.Y()/eofs.Y();

        return hPolynomialThermo<EquationOfState, PolySize>
        (
            eofs,
            Y1*pt1.hf_ + Y2*pt2.hf_,
            Y1*pt1.sf_ + Y2*pt2.sf_,
            Y1*pt1.CpCoeffs_ + Y2*pt2.CpCoeffs_,
            Y1*pt1.hCoeffs_ + Y2*pt2.hCoeffs_,
            Y1*pt1.sCoeffs_ + Y2*pt2.sCoeffs_
        );
    }
}


template<class EquationOfState, int PolySize>
inline Foam::hPolynomialThermo<EquationOfState, PolySize> Foam::operator*
(
    const scalar s,
    const hPolynomialThermo<EquationOfState, PolySize>& pt
)
{
    return hPolynomialThermo<EquationOfState, PolySize>
    (
        s*static_cast<const EquationOfState&>(pt),
        pt.hf_,
        pt.sf_,
        pt.CpCoeffs_,
        pt.hCoeffs_,
        pt.sCoeffs_
    );
}


template<class EquationOfState, int PolySize>
inline Foam::hPolynomialThermo<EquationOfState, PolySize> Foam::operator==
(
    const hPolynomialThermo<EquationOfState, PolySize>& pt1,
    const hPolynomialThermo<EquationOfState, PolySize>& pt2
)
{
    EquationOfState eofs
    (
        static_cast<const EquationOfState&>(pt1)
     == static_cast<const EquationOfState&>(pt2)
    );

    const scalar Y1 = pt1.Y()/eofs.Y();
    const scalar Y2 = pt2.Y()/eofs.Y();

    return hPolynomialThermo<EquationOfState, PolySize>
    (
        eofs,
        Y2*pt2.hf_       - Y1*pt1.hf_,
        Y2*pt2.sf_       - Y1*pt1.sf_,
        Y2*pt2.CpCoeffs_ - Y1*pt1.CpCoeffs_,
        Y2*pt2.hCoeffs_  - Y1*pt1.hCoeffs_,
        Y2*pt2.sCoeffs_  - Y1*pt1.sCoeffs_
    );
}


// ************************************************************************* //
