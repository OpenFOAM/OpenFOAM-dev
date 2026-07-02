/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2020-2026 OpenFOAM Foundation
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

#include "ePolynomialThermo.H"
#include "thermodynamicConstants.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class EquationOfState, int PolySize>
inline Foam::ePolynomialThermo<EquationOfState, PolySize>::ePolynomialThermo
(
    const EquationOfState& pt,
    const scalar hf,
    const scalar sf,
    const FixedLaurentPolynomial<scalar, 0, PolySize>& CvCoeffs,
    const scalar esRef,
    const scalar sRef
)
:
    EquationOfState(pt),
    hf_(hf),
    sf_(sf),
    CvCoeffs_(CvCoeffs),
    esRef_(esRef),
    sRef_(sRef)
{}


template<class EquationOfState, int PolySize>
inline Foam::ePolynomialThermo<EquationOfState, PolySize>::ePolynomialThermo
(
    const word& name,
    const ePolynomialThermo& pt
)
:
    EquationOfState(name, pt),
    hf_(pt.hf_),
    sf_(pt.sf_),
    CvCoeffs_(pt.CvCoeffs_),
    esRef_(pt.esRef_),
    sRef_(pt.sRef_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class EquationOfState, int PolySize>
inline Foam::scalar Foam::ePolynomialThermo<EquationOfState, PolySize>::limit
(
    const scalar T
) const
{
    return T;
}


template<class EquationOfState, int PolySize>
inline Foam::scalar Foam::ePolynomialThermo<EquationOfState, PolySize>::Cv
(
    const scalar p,
    const scalar T
) const
{
    return CvCoeffs_.value(T) + EquationOfState::Cv(p, T);
}


template<class EquationOfState, int PolySize>
inline Foam::scalar Foam::ePolynomialThermo<EquationOfState, PolySize>::es
(
    const scalar p,
    const scalar T
) const
{
    return CvCoeffs_.integral(T) - esRef_ + EquationOfState::e(p, T);
}


template<class EquationOfState, int PolySize>
inline Foam::scalar Foam::ePolynomialThermo<EquationOfState, PolySize>::ea
(
    const scalar p,
    const scalar T
) const
{
    return es(p, T) + hf();
}


template<class EquationOfState, int PolySize>
inline Foam::scalar Foam::ePolynomialThermo<EquationOfState, PolySize>::hf()
const
{
    return hf_;
}


template<class EquationOfState, int PolySize>
inline Foam::scalar Foam::ePolynomialThermo<EquationOfState, PolySize>::s
(
    const scalar p,
    const scalar T
) const
{
    return
        CvCoeffs_.byX().integral(T) - sRef_ + sf_
      + EquationOfState::sv(p, T);
}


template<class EquationOfState, int PolySize>
inline Foam::scalar Foam::ePolynomialThermo<EquationOfState, PolySize>::gStd
(
    const scalar T
) const
{
    using namespace constant::thermodynamic;

    return
        CvCoeffs_.integral(T) - esRef_ + hf()
      + pStd/EquationOfState::rho(pStd, T)
      - s(pStd, T)*T;
}


template<class EquationOfState, int PolySize>
inline Foam::scalar Foam::ePolynomialThermo<EquationOfState, PolySize>::dCpdT
(
    const scalar p,
    const scalar T
) const
{
    NotImplemented;
    return CvCoeffs_.derivative(T); // + EquationOfState::dCpdT
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class EquationOfState, int PolySize>
inline void Foam::ePolynomialThermo<EquationOfState, PolySize>::operator+=
(
    const ePolynomialThermo<EquationOfState, PolySize>& pt
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
        CvCoeffs_ = Y1*CvCoeffs_ + Y2*pt.CvCoeffs_;
        esRef_ = Y1*esRef_ + Y2*pt.esRef_;
        sRef_ = Y1*sRef_ + Y2*pt.sRef_;
    }
}


template<class EquationOfState, int PolySize>
inline void Foam::ePolynomialThermo<EquationOfState, PolySize>::operator*=
(
    const scalar s
)
{
    EquationOfState::operator*=(s);
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template<class EquationOfState, int PolySize>
inline Foam::ePolynomialThermo<EquationOfState, PolySize> Foam::operator+
(
    const ePolynomialThermo<EquationOfState, PolySize>& pt1,
    const ePolynomialThermo<EquationOfState, PolySize>& pt2
)
{
    EquationOfState eofs = pt1;
    eofs += pt2;

    if (mag(eofs.Y()) < small)
    {
        return ePolynomialThermo<EquationOfState, PolySize>
        (
            eofs,
            pt1.hf_,
            pt1.sf_,
            pt1.CvCoeffs_,
            pt1.esRef_,
            pt1.sRef_
        );
    }
    {
        const scalar Y1 = pt1.Y()/eofs.Y();
        const scalar Y2 = pt2.Y()/eofs.Y();

        return ePolynomialThermo<EquationOfState, PolySize>
        (
            eofs,
            Y1*pt1.hf_ + Y2*pt2.hf_,
            Y1*pt1.sf_ + Y2*pt2.sf_,
            Y1*pt1.CvCoeffs_ + Y2*pt2.CvCoeffs_,
            Y1*pt1.esRef_ + Y2*pt2.esRef_,
            Y1*pt1.sRef_ + Y2*pt2.sRef_
        );
    }
}


template<class EquationOfState, int PolySize>
inline Foam::ePolynomialThermo<EquationOfState, PolySize> Foam::operator*
(
    const scalar s,
    const ePolynomialThermo<EquationOfState, PolySize>& pt
)
{
    return ePolynomialThermo<EquationOfState, PolySize>
    (
        s*static_cast<const EquationOfState&>(pt),
        pt.hf_,
        pt.sf_,
        pt.CvCoeffs_,
        pt.esRef_,
        pt.sRef_
    );
}


template<class EquationOfState, int PolySize>
inline Foam::ePolynomialThermo<EquationOfState, PolySize> Foam::operator==
(
    const ePolynomialThermo<EquationOfState, PolySize>& pt1,
    const ePolynomialThermo<EquationOfState, PolySize>& pt2
)
{
    EquationOfState eofs
    (
        static_cast<const EquationOfState&>(pt1)
     == static_cast<const EquationOfState&>(pt2)
    );

    const scalar Y1 = pt1.Y()/eofs.Y();
    const scalar Y2 = pt2.Y()/eofs.Y();

    return ePolynomialThermo<EquationOfState, PolySize>
    (
        eofs,
        Y2*pt2.hf_       - Y1*pt1.hf_,
        Y2*pt2.sf_       - Y1*pt1.sf_,
        Y2*pt2.CvCoeffs_ - Y1*pt1.CvCoeffs_,
        Y2*pt2.esRef_    - Y1*pt1.esRef_,
        Y2*pt2.sRef_     - Y1*pt1.sRef_
    );
}


// ************************************************************************* //
