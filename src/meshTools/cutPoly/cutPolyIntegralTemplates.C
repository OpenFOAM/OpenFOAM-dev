/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2018-2022 OpenFOAM Foundation
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

#include "cutPolyIntegral.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace cutPoly
{

template<class Type>
Type average(const UIndirectList<Type>& l)
{
    Type result = pTraits<Type>::zero;
    forAll(l, i)
    {
        result += l[i];
    }
    return result/l.size();
}

template<class Type>
Type average(const List<Type>& xs, const labelUList& is)
{
    return average(UIndirectList<Type>(xs, is));
}

template<class Container>
typename Container::value_type iterableAverage(const Container& xs)
{
    label n = 0;
    typename Container::value_type nResult =
        pTraits<typename Container::value_type>::zero;

    forAllConstIter(typename Container, xs, iter)
    {
        ++ n;
        nResult += *iter;
    }

    return nResult/n;
}

} // End namespace cutPoly
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace cutPoly
{

template<class Op, class Tuple, class Int, Int ... Is>
auto tupleOp
(
    const Tuple& tuple,
    const std::integer_sequence<Int, Is ...>&,
    const Op& op
)
{
    return std::make_tuple(op(std::get<Is>(tuple)) ...);
}


template<class Op, class ... Types>
auto tupleOp(const std::tuple<Types...>& tuple, const Op& op)
{
    return tupleOp(tuple, std::make_index_sequence<sizeof ... (Types)>(), op);
}


struct OpBegin
{
    template<class Type>
    auto operator()(const Type& x) const
    {
        return x.begin();
    }
};


struct OpDereference
{
    template<class Type>
    auto operator()(const Type& x) const
    {
        return *x;
    }
};


struct OpNext
{
    template<class Type>
    auto operator()(const Type& x) const
    {
        return x.next();
    }
};


template<class ScaleType>
struct OpScaled
{
    const ScaleType s_;

    OpScaled(const ScaleType& s)
    :
        s_(s)
    {}

    template<class Type>
    auto operator()(const Type& x) const
    {
        return s_*x;
    }
};


template<class Op, class Tuple, class Int, Int ... Is>
void tupleInPlaceOp
(
    Tuple& tuple,
    const std::integer_sequence<Int, Is ...>&,
    const Op& op
)
{
    (void)std::initializer_list<nil>
    {(
        op(std::get<Is>(tuple)),
        nil()
    ) ... };
}


template<class Op, class ... Types>
void tupleInPlaceOp(std::tuple<Types...>& tuple, const Op& op)
{
    tupleInPlaceOp(tuple, std::make_index_sequence<sizeof ... (Types)>(), op);
}


struct InPlaceOpAdvance
{
    template<class Type>
    void operator()(Type& x) const
    {
        ++ x;
    }
};


template<class BinaryOp, class Tuple, class Int, Int ... Is>
auto tupleBinaryOp
(
    const Tuple& tupleA,
    const Tuple& tupleB,
    const std::integer_sequence<Int, Is ...>&,
    const BinaryOp& bop
)
{
    return std::make_tuple(bop(std::get<Is>(tupleA), std::get<Is>(tupleB)) ...);
}


template<class BinaryOp, class ... TypesA, class ... TypesB>
auto tupleBinaryOp
(
    const std::tuple<TypesA ...>& tupleA,
    const std::tuple<TypesB ...>& tupleB,
    const BinaryOp& bop
)
{
    return tupleBinaryOp
    (
        tupleA,
        tupleB,
        std::make_index_sequence<sizeof ... (TypesA)>(),
        bop
    );
}


struct BinaryOpAdd
{
    template<class Type>
    auto operator()(const Type& a, const Type& b) const
    {
        return a + b;
    }
};

} // End namespace cutPoly
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<template<class> class FaceValues, class ... Types>
std::tuple
<
    Foam::vector,
    typename Foam::outerProduct<Foam::vector, Types>::type ...
>
Foam::cutPoly::faceAreaIntegral
(
    const FaceValues<point>& fPs,
    const point& fPAvg,
    const std::tuple<FaceValues<Types> ...>& fPsis,
    const std::tuple<Types ...>& fPsiAvg
)
{
    vector fCutsArea = Zero;
    auto fCutsAreaPsi =
        std::make_tuple(typename outerProduct<vector, Types>::type(Zero) ...);

    typename FaceValues<point>::const_iterator fPIter(fPs.begin());
    auto fPsiIter = tupleOp(fPsis, OpBegin());

    for(; fPIter != fPs.end(); ++ fPIter)
    {
        const point p0 = *fPIter;
        const point p1 = fPIter.next();
        auto psi0 = tupleOp(fPsiIter, OpDereference());
        auto psi1 = tupleOp(fPsiIter, OpNext());

        const vector a = ((p1 - p0)^(fPAvg - p0))/2;
        auto v =
            tupleBinaryOp
            (
                psi0,
                tupleBinaryOp
                (
                    psi1,
                    fPsiAvg,
                    BinaryOpAdd()
                ),
                BinaryOpAdd()
            );

        auto av = tupleOp(v, OpScaled<vector>(a/3));

        fCutsArea += a;
        fCutsAreaPsi = tupleBinaryOp(fCutsAreaPsi, av, BinaryOpAdd());

        tupleInPlaceOp(fPsiIter, InPlaceOpAdvance());
    }

    return std::tuple_cat(std::tuple<vector>(fCutsArea), fCutsAreaPsi);
}


template<class Type>
Foam::Tuple2
<
    Foam::vector,
    typename Foam::outerProduct<Foam::vector, Type>::type
>
Foam::cutPoly::faceCutAreaIntegral
(
    const face& f,
    const vector& fArea,
    const Type& fPsi,
    const List<labelPair>& fCuts,
    const pointField& ps,
    const Field<Type>& pPsis,
    const scalarField& pAlphas,
    const scalar isoAlpha,
    const bool below
)
{
    typedef typename outerProduct<vector, Type>::type IntegralType;

    // If there are no cuts return either the entire face or zero, depending on
    // which side of the iso-surface the face is
    if (fCuts.size() == 0)
    {
        if ((pAlphas[f[0]] < isoAlpha) == below)
        {
            return Tuple2<vector, IntegralType>(fArea, fArea*fPsi);
        }
        else
        {
            return Tuple2<vector, IntegralType>(Zero, Zero);
        }
    }

    auto result =
        faceAreaIntegral
        (
            cutPoly::FaceCutValues<vector>
            (
                f,
                fCuts,
                ps,
                pAlphas,
                isoAlpha,
                below
            ),
            average(ps, f),
            std::make_tuple
            (
                cutPoly::FaceCutValues<Type>
                (
                    f,
                    fCuts,
                    pPsis,
                    pAlphas,
                    isoAlpha,
                    below
                )
            ),
            std::make_tuple(average(pPsis, f))
        );

    return
        Tuple2<vector, IntegralType>
        (
            std::get<0>(result),
            std::get<1>(result)
        );
}


template<class Type>
Foam::Tuple2<Foam::scalar, Type> Foam::cutPoly::cellCutVolumeIntegral
(
    const cell& c,
    const cellEdgeAddressing& cAddr,
    const scalar cVolume,
    const Type& cPsi,
    const labelListList& cCuts,
    const faceUList& fs,
    const vectorField& fAreas,
    const vectorField& fCentres,
    const vectorField& fPsis,
    const vectorField& fCutAreas,
    const vectorField& fCutPsis,
    const pointField& ps,
    const Field<Type>& pPsis,
    const scalarField& pAlphas,
    const scalar isoAlpha,
    const bool below
)
{
    // If there are no cuts return either the entire cell or zero, depending on
    // which side of the iso-surface the cell is
    if (cCuts.size() == 0)
    {
        if ((pAlphas[fs[c[0]][0]] < isoAlpha) == below)
        {
            return Tuple2<scalar, Type>(cVolume, cVolume*cPsi);
        }
        else
        {
            return Tuple2<scalar, Type>(0, Zero);
        }
    }

    // Averages
    const point cPAvg = average(fCentres, c);
    const Type cPsiAvg = average(fPsis, c);

    // Initialise totals
    scalar cCutsVolume = 0;
    Type cCutsVolumePsi = Zero;

    // Face contributions
    forAll(c, cfi)
    {
        // We use the un-cut face's centroid as the base of the pyramid formed
        // by the cut face. This is potentially less exact, but it is
        // consistent with the un-cut cell volume calculation. This means if
        // you run this function with both values of "below" then the result
        // will exactly sum to the volume of the overall cell. If we used the
        // cut-face's centroid this would not be the case.
        const vector& fBaseP = fCentres[c[cfi]];

        const scalar pyrVolume =
            (cAddr.cOwns()[cfi] ? +1 : -1)
           *(fCutAreas[c[cfi]] & (fBaseP - cPAvg))/3;

        cCutsVolume += pyrVolume;
        cCutsVolumePsi += pyrVolume*(3*fCutPsis[c[cfi]] + cPsiAvg)/4;
    }

    // Cut contributions
    {
        const cutPoly::CellCutValues<point> cPValues
        (
            c,
            cAddr,
            cCuts,
            fs,
            ps,
            pAlphas,
            isoAlpha
        );

        const cutPoly::CellCutValues<Type> cPsiValues
        (
            c,
            cAddr,
            cCuts,
            fs,
            pPsis,
            pAlphas,
            isoAlpha
        );

        // This method is more exact, as it uses the true centroid of the cell
        // cut. However, to obtain that centroid we have to divide by the area
        // magnitude, so this can't be generalised to types that do not support
        // division.
        /*
        auto fSumPPsis =
            faceAreaIntegral
            (
                cPValues,
                iterableAverage(cPValues),
                std::make_tuple
                (
                    cPValues,
                    cPsiValues
                ),
                std::make_tuple
                (
                    iterableAverage(cPValues),
                    iterableAverage(cPsiValues)
                )
            );

        const vector& fCutArea = std::get<0>(fSumPPsis);
        const scalar fMagSqrCutArea = magSqr(fCutArea);

        const point fCutCentre =
            fMagSqrCutArea > vSmall
          ? (fCutArea & std::get<1>(fSumPPsis))/fMagSqrCutArea
          : cPAvg;
        const Type fCutPsi =
            fMagSqrCutArea > vSmall
          ? (fCutArea & std::get<2>(fSumPPsis))/fMagSqrCutArea
          : cPsiAvg;
        */

        // This method is more approximate, as it uses point averages rather
        // than centroids. This does not involve division, though, so this
        // could be used with types like polynomials that only support addition
        // and multiplication.
        auto fSumPPsis =
            faceAreaIntegral
            (
                cPValues,
                iterableAverage(cPValues),
                std::make_tuple(),
                std::make_tuple()
            );

        const vector& fCutArea = std::get<0>(fSumPPsis);
        const point fCutCentre = iterableAverage(cPValues);
        const Type fCutPsi = iterableAverage(cPsiValues);

        const scalar pyrVolume =
            (below ? +1 : -1)*(fCutArea & (fCutCentre - cPAvg))/3;

        cCutsVolume += pyrVolume;
        cCutsVolumePsi += pyrVolume*(3*fCutPsi + cPsiAvg)/4;
    }

    return Tuple2<scalar, Type>(cCutsVolume, cCutsVolumePsi);
}


// ************************************************************************* //
