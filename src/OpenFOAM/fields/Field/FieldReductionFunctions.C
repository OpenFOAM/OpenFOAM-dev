/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2026 OpenFOAM Foundation
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

#include "PstreamReduceOps.H"
#include "FieldM.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define G_REDUCTION_FUNCTION(ReturnType, gFunc, Func, rFunc)                   \
                                                                               \
    template<class Type>                                                       \
    ReturnType gFunc                                                           \
    (                                                                          \
        const UList<Type>& f,                                                  \
        const label comm                                                       \
    )                                                                          \
    {                                                                          \
        ReturnType res = Func(f);                                              \
        reduce(res, rFunc##Op(), Pstream::msgType(), comm);              \
        return res;                                                            \
    }

#define TMP_REDUCTION_FUNCTION(ReturnType, Func)                               \
                                                                               \
    template<class Type>                                                       \
    ReturnType Func(const tmp<Field<Type>>& tf1)                               \
    {                                                                          \
        ReturnType res = Func(tf1());                                          \
        tf1.clear();                                                           \
        return res;                                                            \
    }

/* * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * */

template<class Type>
Type max(const UList<Type>& f)
{
    if (f.size())
    {
        Type Max(f[0]);
        TFOR_ALL_S_OP_FUNC_F_S(Type, Max, =, max, Type, f, Type, Max)
        return Max;
    }
    else
    {
        return pTraits<Type>::min;
    }
}

G_REDUCTION_FUNCTION(Type, gMax, max, max)
TMP_REDUCTION_FUNCTION(Type, max)
TMP_REDUCTION_FUNCTION(Type, gMax)

template<class Type>
Type min(const UList<Type>& f)
{
    if (f.size())
    {
        Type Min(f[0]);
        TFOR_ALL_S_OP_FUNC_F_S(Type, Min, =, min, Type, f, Type, Min)
        return Min;
    }
    else
    {
        return pTraits<Type>::max;
    }
}

G_REDUCTION_FUNCTION(Type, gMin, min, min)
TMP_REDUCTION_FUNCTION(Type, min)
TMP_REDUCTION_FUNCTION(Type, gMin)

template<class Type>
Type sum(const UList<Type>& f)
{
    if (f.size())
    {
        Type Sum = Zero;
        TFOR_ALL_S_OP_F(Type, Sum, +=, Type, f)
        return Sum;
    }
    else
    {
        return Zero;
    }
}

G_REDUCTION_FUNCTION(Type, gSum, sum, sum)
TMP_REDUCTION_FUNCTION(Type, sum)
TMP_REDUCTION_FUNCTION(Type, gSum)

template<class Type>
Type maxMagSqr(const UList<Type>& f)
{
    if (f.size())
    {
        Type Max(f[0]);
        TFOR_ALL_S_OP_FUNC_F_S
        (
            Type,
            Max,
            =,
            maxMagSqrOp(),
            Type,
            f,
            Type,
            Max
        )
        return Max;
    }
    else
    {
        return Zero;
    }
}

G_REDUCTION_FUNCTION(Type, gMaxMagSqr, maxMagSqr, maxMagSqr)
TMP_REDUCTION_FUNCTION(Type, maxMagSqr)
TMP_REDUCTION_FUNCTION(Type, gMaxMagSqr)

template<class Type>
Type minMagSqr(const UList<Type>& f)
{
    if (f.size())
    {
        Type Min(f[0]);
        TFOR_ALL_S_OP_FUNC_F_S
        (
            Type,
            Min,
            =,
            minMagSqrOp(),
            Type,
            f,
            Type,
            Min
        )
        return Min;
    }
    else
    {
        return pTraits<Type>::rootMax;
    }
}

G_REDUCTION_FUNCTION(Type, gMinMagSqr, minMagSqr, minMagSqr)
TMP_REDUCTION_FUNCTION(Type, minMagSqr)
TMP_REDUCTION_FUNCTION(Type, gMinMagSqr)

template<class Type>
scalar sumProd(const UList<Type>& f1, const UList<Type>& f2)
{
    if (f1.size() && (f1.size() == f2.size()))
    {
        scalar SumProd = 0;
        TFOR_ALL_S_OP_F_OP_F(scalar, SumProd, +=, Type, f1, &&, Type, f2)
        return SumProd;
    }
    else
    {
        return 0;
    }
}

template<class Type>
scalar gSumProd
(
    const UList<Type>& f1,
    const UList<Type>& f2,
    const label comm
)
{
    scalar SumProd = sumProd(f1, f2);
    reduce(SumProd, sumOp(), Pstream::msgType(), comm);
    return SumProd;
}

template<class Type>
Type sumCmptProd(const UList<Type>& f1, const UList<Type>& f2)
{
    if (f1.size() && (f1.size() == f2.size()))
    {
        Type SumProd = Zero;
        TFOR_ALL_S_OP_FUNC_F_F
        (
            Type,
            SumProd,
            +=,
            cmptMultiply,
            Type,
            f1,
            Type,
            f2
        )
        return SumProd;
    }
    else
    {
        return Zero;
    }
}

template<class Type>
Type gSumCmptProd
(
    const UList<Type>& f1,
    const UList<Type>& f2,
    const label comm
)
{
    Type SumProd = sumCmptProd(f1, f2);
    reduce(SumProd, sumOp(), Pstream::msgType(), comm);
    return SumProd;
}

template<class Type>
scalar sumSqr(const UList<Type>& f)
{
    if (f.size())
    {
        scalar SumSqr = 0;
        TFOR_ALL_S_OP_FUNC_F(scalar, SumSqr, +=, sqr, Type, f)
        return SumSqr;
    }
    else
    {
        return 0;
    }
}

G_REDUCTION_FUNCTION(scalar, gSumSqr, sumSqr, sum)
TMP_REDUCTION_FUNCTION(scalar, sumSqr)
TMP_REDUCTION_FUNCTION(scalar, gSumSqr)

template<class Type>
scalar sumMag(const UList<Type>& f)
{
    if (f.size())
    {
        scalar SumMag = 0;
        TFOR_ALL_S_OP_FUNC_F(scalar, SumMag, +=, mag, Type, f)
        return SumMag;
    }
    else
    {
        return 0;
    }
}

G_REDUCTION_FUNCTION(scalar, gSumMag, sumMag, sum)
TMP_REDUCTION_FUNCTION(scalar, sumMag)
TMP_REDUCTION_FUNCTION(scalar, gSumMag)

template<class Type>
Type sumCmptMag(const UList<Type>& f)
{
    if (f.size())
    {
        Type SumMag = Zero;
        TFOR_ALL_S_OP_FUNC_F(scalar, SumMag, +=, cmptMag, Type, f)
        return SumMag;
    }
    else
    {
        return Zero;
    }
}

G_REDUCTION_FUNCTION(Type, gSumCmptMag, sumCmptMag, sum)
TMP_REDUCTION_FUNCTION(Type, sumCmptMag)
TMP_REDUCTION_FUNCTION(Type, gSumCmptMag)

template<class Type>
Type average(const UList<Type>& f)
{
    if (f.size())
    {
        Type avrg = sum(f)/f.size();

        return avrg;
    }
    else
    {
        WarningInFunction
            << "empty field, returning zero" << endl;

        return Zero;
    }
}

template<class Type>
Type gAverage
(
    const UList<Type>& f,
    const label comm
)
{
    label n = f.size();
    Type s = sum(f);
    sumReduce(s, n, Pstream::msgType(), comm);

    if (n > 0)
    {
        Type avrg = s/n;

        return avrg;
    }
    else
    {
        WarningInFunction
            << "empty field, returning zero." << endl;

        return Zero;
    }
}

TMP_REDUCTION_FUNCTION(Type, average)
TMP_REDUCTION_FUNCTION(Type, gAverage)

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#undef G_REDUCTION_FUNCTION
#undef TMP_REDUCTION_FUNCTION

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
