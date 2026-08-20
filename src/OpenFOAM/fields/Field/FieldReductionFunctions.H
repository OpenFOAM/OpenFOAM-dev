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

#include "UPstream.H"

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
        const label comm = UPstream::worldComm                                 \
    );                                                                         \
                                                                               \
    template<class Type>                                                       \
    ReturnType gFunc                                                           \
    (                                                                          \
        const tmp<Field<Type>>& tf,                                            \
        const label comm                                                       \
    );

#define TMP_REDUCTION_FUNCTION(ReturnType, Func)                               \
                                                                               \
    template<class Type>                                                       \
    ReturnType Func(const tmp<Field<Type>>& tf1);

// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

template<class Type>
Type max(const UList<Type>& f);

G_REDUCTION_FUNCTION(Type, gMax, max, max)
TMP_REDUCTION_FUNCTION(Type, max)
TMP_REDUCTION_FUNCTION(Type, gMax)

template<class Type>
Type min(const UList<Type>& f);

G_REDUCTION_FUNCTION(Type, gMin, min, min)
TMP_REDUCTION_FUNCTION(Type, min)
TMP_REDUCTION_FUNCTION(Type, gMin)

template<class Type>
Type sum(const UList<Type>& f);

G_REDUCTION_FUNCTION(Type, gSum, sum, sum)
TMP_REDUCTION_FUNCTION(Type, sum)
TMP_REDUCTION_FUNCTION(Type, gSum)

template<class Type>
Type maxMagSqr(const UList<Type>& f);

G_REDUCTION_FUNCTION(Type, gMaxMagSqr, maxMagSqr, maxMagSqr)
TMP_REDUCTION_FUNCTION(Type, maxMagSqr)
TMP_REDUCTION_FUNCTION(Type, gMaxMagSqr)

template<class Type>
Type minMagSqr(const UList<Type>& f);

G_REDUCTION_FUNCTION(Type, gMinMagSqr, minMagSqr, minMagSqr)
TMP_REDUCTION_FUNCTION(Type, minMagSqr)
TMP_REDUCTION_FUNCTION(Type, gMinMagSqr)

template<class Type>
scalar sumProd(const UList<Type>& f1, const UList<Type>& f2);

template<class Type>
scalar gSumProd
(
    const UList<Type>& f1,
    const UList<Type>& f2,
    const label comm = UPstream::worldComm
);

template<class Type>
Type sumCmptProd(const UList<Type>& f1, const UList<Type>& f2);

template<class Type>
Type gSumCmptProd
(
    const UList<Type>& f1,
    const UList<Type>& f2,
    const label comm = UPstream::worldComm
);

template<class Type>
scalar sumSqr(const UList<Type>& f);

G_REDUCTION_FUNCTION(scalar, gSumSqr, sumSqr, sum)
TMP_REDUCTION_FUNCTION(scalar, sumSqr)
TMP_REDUCTION_FUNCTION(scalar, gSumSqr)

template<class Type>
scalar sumMag(const UList<Type>& f);

G_REDUCTION_FUNCTION(scalar, gSumMag, sumMag, sum)
TMP_REDUCTION_FUNCTION(scalar, sumMag)
TMP_REDUCTION_FUNCTION(scalar, gSumMag)

template<class Type>
Type sumCmptMag(const UList<Type>& f);

G_REDUCTION_FUNCTION(Type, gSumCmptMag, sumCmptMag, sum)
TMP_REDUCTION_FUNCTION(Type, sumCmptMag)
TMP_REDUCTION_FUNCTION(Type, gSumCmptMag)

template<class Type>
Type average(const UList<Type>& f);

TMP_REDUCTION_FUNCTION(Type, average)

template<class Type>
Type gAverage
(
    const UList<Type>& f,
    const label comm = UPstream::worldComm
);

TMP_REDUCTION_FUNCTION(Type, gAverage)

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#undef G_REDUCTION_FUNCTION
#undef TMP_REDUCTION_FUNCTION

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
