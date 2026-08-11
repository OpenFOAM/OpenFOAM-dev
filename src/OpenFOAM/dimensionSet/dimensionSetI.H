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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

inline scalar dimensionSet::operator[](const dimensionType type) const
{
    return exponents_[type];
}


inline scalar& dimensionSet::operator[](const dimensionType type)
{
    return exponents_[type];
}


inline scalar dimensionSet::operator[](const label type) const
{
    return exponents_[type];
}


inline scalar& dimensionSet::operator[](const label type)
{
    return exponents_[type];
}


inline bool dimensionSet::operator!=(const dimensionSet& ds) const
{
    return !(operator==(ds));
}


inline bool dimensionSet::operator=(const dimensionSet& ds) const
{
    if (dimensionSet::debug && *this != ds)
    {
        FatalErrorInFunction
            << "Different dimensions for =" << endl
            << "     dimensions : " << info() << " = " << ds.info() << endl
            << abort(FatalError);
    }

    return true;
}


inline bool dimensionSet::operator+=(const dimensionSet& ds) const
{
    if (dimensionSet::debug && *this != ds)
    {
        FatalErrorInFunction
            << "Different dimensions for +=" << endl
            << "     dimensions : " << info() << " = " << ds.info() << endl
            << abort(FatalError);
    }

    return true;
}


inline bool dimensionSet::operator-=(const dimensionSet& ds) const
{
    if (dimensionSet::debug && *this != ds)
    {
        FatalErrorInFunction
            << "Different dimensions for -=" << endl
            << "     dimensions : " << info() << " = " << ds.info() << endl
            << abort(FatalError);
    }

    return true;
}


inline bool dimensionSet::operator*=(const dimensionSet& ds)
{
    reset((*this)*ds);

    return true;
}


inline bool dimensionSet::operator/=(const dimensionSet& ds)
{
    reset((*this)/ds);

    return true;
}


// * * * * * * * * * * * * * * * Friend functions * * * * * * * * * * * * * * //

inline dimensionSet max(const dimensionSet& ds1, const dimensionSet& ds2)
{
    if (dimensionSet::debug && ds1 != ds2)
    {
        FatalErrorInFunction
            << "Arguments of max have different dimensions" << endl
            << "     dimensions : " << ds1.info() << " and " << ds2.info()
            << endl << abort(FatalError);
    }

    return ds1;
}


inline dimensionSet min(const dimensionSet& ds1, const dimensionSet& ds2)
{
    if (dimensionSet::debug && ds1 != ds2)
    {
        FatalErrorInFunction
            << "Arguments of min have different dimensions" << endl
            << "     dimensions : " << ds1.info() << " and " << ds2.info()
            << endl << abort(FatalError);
    }

    return ds1;
}


inline dimensionSet cmptMultiply
(
    const dimensionSet& ds1,
    const dimensionSet& ds2
)
{
    return ds1*ds2;
}


inline dimensionSet cmptDivide
(
    const dimensionSet& ds1,
    const dimensionSet& ds2
)
{
    return ds1/ds2;
}


inline dimensionSet cmptMag(const dimensionSet& ds)
{
    return ds;
}


inline dimensionSet pow(const dimensionSet& ds, const scalar p)
{
    dimensionSet dimPow
    (
        ds[dimensionSet::MASS]*p,
        ds[dimensionSet::LENGTH]*p,
        ds[dimensionSet::TIME]*p,
        ds[dimensionSet::TEMPERATURE]*p,
        ds[dimensionSet::MOLES]*p,
        ds[dimensionSet::CURRENT]*p,
        ds[dimensionSet::LUMINOUS_INTENSITY]*p
    );

    return dimPow;
}


inline dimensionSet sqr(const dimensionSet& ds)
{
    return pow(ds, 2);
}


inline dimensionSet innerSqr(const dimensionSet& ds)
{
    return sqr(ds);
}


inline dimensionSet pow3(const dimensionSet& ds)
{
    return pow(ds, 3);
}


inline dimensionSet pow4(const dimensionSet& ds)
{
    return pow(ds, 4);
}


inline dimensionSet pow5(const dimensionSet& ds)
{
    return pow(ds, 5);
}


inline dimensionSet pow6(const dimensionSet& ds)
{
    return pow(ds, 6);
}


inline dimensionSet sqrt(const dimensionSet& ds)
{
    return pow(ds, 0.5);
}


inline dimensionSet pow025(const dimensionSet& ds)
{
    return sqrt(sqrt(ds));
}


inline dimensionSet cbrt(const dimensionSet& ds)
{
    return pow(ds, 1.0/3.0);
}


inline dimensionSet magSqr(const dimensionSet& ds)
{
    return pow(ds, 2);
}


inline dimensionSet mag(const dimensionSet& ds)
{
    return ds;
}


inline dimensionSet sign(const dimensionSet&)
{
    return dimless;
}


inline dimensionSet pos(const dimensionSet&)
{
    return dimless;
}


inline dimensionSet pos0(const dimensionSet&)
{
    return dimless;
}


inline dimensionSet neg(const dimensionSet&)
{
    return dimless;
}


inline dimensionSet neg0(const dimensionSet&)
{
    return dimless;
}


inline dimensionSet posPart(const dimensionSet& ds)
{
    return ds;
}


inline dimensionSet negPart(const dimensionSet& ds)
{
    return ds;
}


#define transFunc(Func)                                                        \
    inline dimensionSet Func(const dimensionSet& ds)                           \
    {                                                                          \
        if (dimensionSet::debug && !ds.dimensionless())                        \
        {                                                                      \
            FatalErrorInFunction                                               \
                << "Argument of trancendental function " #Func                 \
                   " is not dimensionless"                                     \
                << abort(FatalError);                                          \
        }                                                                      \
                                                                               \
        return ds;                                                             \
    }                                                                          \

transFunc(exp)
transFunc(log)
transFunc(log10)
transFunc(sin)
transFunc(cos)
transFunc(tan)
transFunc(asin)
transFunc(acos)
transFunc(atan)
transFunc(sinh)
transFunc(cosh)
transFunc(tanh)
transFunc(asinh)
transFunc(acosh)
transFunc(atanh)
transFunc(erf)
transFunc(erfc)
transFunc(lgamma)
transFunc(j0)
transFunc(j1)
transFunc(y0)
transFunc(y1)

#undef transFunc


inline dimensionSet atan2(const dimensionSet& ds1, const dimensionSet& ds2)
{
    if (dimensionSet::debug && ds1 != ds2)
    {
        FatalErrorInFunction
            << "Arguments of atan2 have different dimensions" << endl
            << "     dimensions : " << ds1.info() << " and " << ds2.info()
            << endl << abort(FatalError);
    }

    return dimless;
}


inline dimensionSet jn(const dimensionSet& dn, const dimensionSet& ds)
{
    if (dimensionSet::debug && !dn.dimensionless() && !ds.dimensionless())
    {
        FatalErrorInFunction
            << "Arguments of jn are not dimensionless" << endl
            << "     dimensions : " << dn.info() << " and " << ds.info()
            << endl << abort(FatalError);
    }

    return dimless;
}


inline dimensionSet yn(const dimensionSet& dn, const dimensionSet& ds)
{
    if (dimensionSet::debug && !dn.dimensionless() && !ds.dimensionless())
    {
        FatalErrorInFunction
            << "Arguments of yn are not dimensionless" << endl
            << "     dimensions : " << dn.info() << " and " << ds.info()
            << endl << abort(FatalError);
    }

    return dimless;
}


inline dimensionSet tr(const dimensionSet& ds)
{
    return ds;
}


inline dimensionSet sph(const dimensionSet& ds)
{
    return ds;
}


inline dimensionSet diag(const dimensionSet& ds)
{
    return ds;
}


inline dimensionSet symm(const dimensionSet& ds)
{
    return ds;
}


inline dimensionSet twoSymm(const dimensionSet& ds)
{
    return ds;
}


inline dimensionSet dev(const dimensionSet& ds)
{
    return ds;
}


inline dimensionSet dev2(const dimensionSet& ds)
{
    return ds;
}


inline dimensionSet skew(const dimensionSet& ds)
{
    return ds;
}


inline dimensionSet hdual(const dimensionSet& ds)
{
    return ds;
}


inline dimensionSet eigenValues(const dimensionSet& ds)
{
    return ds;
}


inline dimensionSet eigenVectors(const dimensionSet& ds)
{
    return dimless;
}


inline dimensionSet det(const dimensionSet& ds)
{
    return pow3(ds);
}


inline dimensionSet cof(const dimensionSet& ds)
{
    return sqr(ds);
}


inline dimensionSet inv(const dimensionSet& ds)
{
    return dimless/ds;
}


inline dimensionSet normalised(const dimensionSet& ds)
{
    return dimless;
}


inline dimensionSet perpendicular(const dimensionSet& ds)
{
    return dimless;
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

inline dimensionSet operator-(const dimensionSet& ds)
{
    return ds;
}


inline dimensionSet operator+
(
    const dimensionSet& ds1,
    const dimensionSet& ds2
)
{
    dimensionSet dimSum(ds1);

    if (dimensionSet::debug && ds1 != ds2)
    {
        FatalErrorInFunction
            << "LHS and RHS of + have different dimensions" << endl
            << "     dimensions : " << ds1.info() << " + " << ds2.info()
            << endl << abort(FatalError);
    }

    return dimSum;
}


inline dimensionSet operator-
(
    const dimensionSet& ds1,
    const dimensionSet& ds2
)
{
    dimensionSet dimDifference(ds1);

    if (dimensionSet::debug && ds1 != ds2)
    {
        FatalErrorInFunction
            << "LHS and RHS of - have different dimensions" << endl
            << "     dimensions : " << ds1.info() << " - " << ds2.info()
            << endl << abort(FatalError);
    }

    return dimDifference;
}


inline dimensionSet operator&
(
    const dimensionSet& ds1,
    const dimensionSet& ds2
)
{
    return ds1*ds2;
}


inline dimensionSet operator^
(
    const dimensionSet& ds1,
    const dimensionSet& ds2
)
{
    return ds1*ds2;
}


inline dimensionSet operator&&
(
    const dimensionSet& ds1,
    const dimensionSet& ds2
)
{
    return ds1*ds2;
}


inline dimensionSet trans(const dimensionSet& ds)
{
    if (dimensionSet::debug && !ds.dimensionless())
    {
        FatalErrorInFunction
            << "Argument of trancendental function not dimensionless"
            << abort(FatalError);
    }

    return ds;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
