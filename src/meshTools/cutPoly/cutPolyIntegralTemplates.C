/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2018-2025 OpenFOAM Foundation
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
#include <utility>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace cutPoly
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Int, Int Offset, class Sequence>
struct OffsetSequence;

template<class Int, Int Offset, Int ... Is>
struct OffsetSequence<Int, Offset, std::integer_sequence<Int, Is ...>>
{
    using type = std::integer_sequence<Int, Is + Offset ...>;
};


template<class Int, Int Min, Int Max>
struct RangeSequence
{
    using type =
        typename OffsetSequence
        <
            Int,
            Min,
            std::make_index_sequence<Max - Min>
        >::type;
};


template<class Tuple, class Int, Int ... Is>
auto tupleSubset
(
    const Tuple& tuple,
    const std::integer_sequence<Int, Is ...>&
)
{
    return std::make_tuple(std::get<Is>(tuple) ...);
}


template<class Op, class Tuple, class Int, Int ... Is>
auto tupleSubset
(
    const Tuple& tuple,
    const std::integer_sequence<Int, Is ...>&,
    const Op& op
)
{
    return std::make_tuple(op(std::get<Is>(tuple)) ...);
}


template<class ... Types>
auto tupleTail(const std::tuple<Types ...>& tuple)
{
    return
        tupleSubset
        (
            tuple,
            typename RangeSequence<std::size_t, 1, sizeof ... (Types)>::type()
        );
}


template<class Op, class ... Types>
auto tupleOp(const std::tuple<Types ...>& tuple, const Op& op)
{
    return
        tupleSubset
        (
            tuple,
            std::make_index_sequence<sizeof ... (Types)>(),
            op
        );
}


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
void tupleInPlaceOp(std::tuple<Types ...>& tuple, const Op& op)
{
    tupleInPlaceOp(tuple, std::make_index_sequence<sizeof ... (Types)>(), op);
}


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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace cutPoly
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<template<class> class FaceValues, class ... Types>
Foam::Tuple2
<
    Foam::vector,
    std::tuple<Foam::cutPoly::AreaIntegralType<Types> ...>
>
Foam::cutPoly::faceAreaIntegral
(
    const FaceValues<point>& fPs,
    const point& fPAvg,
    const std::tuple<FaceValues<Types> ...>& fPsis,
    const std::tuple<Types ...>& fPsiAvgs
)
{
    vector fCutsArea = Zero;
    auto fCutsAreaPsi = std::make_tuple(AreaIntegralType<Types>(Zero) ...);

    typename FaceValues<point>::const_iterator fPIter(fPs.begin());
    auto fPsiIter = tupleOp(fPsis, OpBegin());

    for(; fPIter != fPs.end(); ++ fPIter)
    {
        const point p0 = *fPIter;
        const point p1 = fPIter.next();
        auto psi0 = tupleOp(fPsiIter, OpDereference());
        auto psi1 = tupleOp(fPsiIter, OpNext());

        const vector a = ((p1 - p0)^(fPAvg - p0))/2;

        fCutsArea += a;
        fCutsAreaPsi =
            tupleBinaryOp
            (
                fCutsAreaPsi,
                tupleOp
                (
                    tupleBinaryOp
                    (
                        psi0,
                        tupleBinaryOp
                        (
                            psi1,
                            fPsiAvgs,
                            BinaryOpAdd()
                        ),
                        BinaryOpAdd()
                    ),
                    OpScaled<vector>(a/3)
                ),
                BinaryOpAdd()
            );

        tupleInPlaceOp(fPsiIter, InPlaceOpAdvance());
    }

    return
        Tuple2<vector, std::tuple<AreaIntegralType<Types> ...>>
        (
            fCutsArea,
            fCutsAreaPsi
        );
}


template<template<class> class FaceValues, class ... Types>
Foam::Tuple2<Foam::vector, std::tuple<Types ...>>
Foam::cutPoly::faceAreaAverage
(
    const FaceValues<point>& fPs,
    const point& fPAvg,
    const std::tuple<FaceValues<Types> ...>& fPsis,
    const std::tuple<Types ...>& fPsiAvgs
)
{
    auto fSumPPsis = faceAreaIntegral(fPs, fPAvg, fPsis, fPsiAvgs);

    const vector& fArea = fSumPPsis.first();
    const scalar fMagSqrArea = magSqr(fArea);

    return
        Tuple2<vector, std::tuple<Types ...>>
        (
            fArea,
            fMagSqrArea > vSmall
          ? tupleOp(fSumPPsis.second(), OpPreInner(fArea/fMagSqrArea))
          : fPsiAvgs
        );
}


template<class ... Types>
Foam::Tuple2
<
    Foam::vector,
    std::tuple<Foam::cutPoly::AreaIntegralType<Types> ...>
>
Foam::cutPoly::faceCutAreaIntegral
(
    const face& f,
    const vector& fArea,
    const std::tuple<Types ...>& fPsis,
    const List<labelPair>& fCuts,
    const pointField& ps,
    const std::tuple<const Field<Types>& ...>& pPsis,
    const scalarField& pAlphas,
    const scalar isoAlpha,
    const bool below
)
{
    // If there are no cuts return either the entire face or zero, depending on
    // which side of the iso-surface the face is
    if (fCuts.size() == 0)
    {
        if ((pAlphas[f[0]] < isoAlpha) == below)
        {
            return
                Tuple2<vector, std::tuple<AreaIntegralType<Types> ...>>
                (
                    fArea,
                    tupleOp(fPsis, OpScaled<vector>(fArea))
                );
        }
        else
        {
            return
                Tuple2<vector, std::tuple<AreaIntegralType<Types> ...>>
                (
                    vector::zero,
                    std::make_tuple(pTraits<AreaIntegralType<Types>>::zero ...)
                );
        }
    }

    return
        faceAreaIntegral
        (
            FaceCutValues<vector>(f, fCuts, ps, pAlphas, isoAlpha, below),
            OpIndirectAverage(f)(ps),
            tupleOp(pPsis, OpFaceCutValues(f, fCuts, pAlphas, isoAlpha, below)),
            tupleOp(pPsis, OpIndirectAverage(f))
        );
}


template<class ... Types>
Foam::Tuple2<Foam::scalar, std::tuple<Types ...>>
Foam::cutPoly::cellVolumeIntegral
(
    const cell& c,
    const cellEdgeAddressing& cAddr,
    const point& cPAvg,
    const std::tuple<Types ...>& cPsiAvgs,
    const vectorField& fAreas,
    const pointField& fCentres,
    const std::tuple<const Field<Types>& ...>& fPsis
)
{
    scalar cVolume = 0;
    std::tuple<Types ...> cVolumePsis(pTraits<Types>::zero ...);

    forAll(c, cfi)
    {
        const scalar pyrVolume =
            (cAddr.cOwns()[cfi] ? +1 : -1)
           *(fAreas[c[cfi]] & (fCentres[c[cfi]] - cPAvg))/3;

        cVolume += pyrVolume;
        cVolumePsis =
            tupleBinaryOp
            (
                cVolumePsis,
                tupleBinaryOp
                (
                    tupleOp
                    (
                        tupleOp(fPsis, OpIndex(c[cfi])),
                        OpScaled<scalar>(scalar(3)/scalar(4)*pyrVolume)
                    ),
                    tupleOp
                    (
                        cPsiAvgs,
                        OpScaled<scalar>(scalar(1)/scalar(4)*pyrVolume)
                    ),
                    BinaryOpAdd()
                ),
                BinaryOpAdd()
            );
    }

    return Tuple2<scalar, std::tuple<Types ...>>(cVolume, cVolumePsis);
}


template<class ... Types>
Foam::Tuple2<Foam::scalar, std::tuple<Types ...>>
Foam::cutPoly::cellCutVolumeIntegral
(
    const cell& c,
    const cellEdgeAddressing& cAddr,
    const scalar cVolume,
    const std::tuple<Types ...>& cPsis,
    const labelListList& cCuts,
    const faceUList& fs,
    const vectorField& fAreas,
    const pointField& fCentres,
    const std::tuple<const Field<Types>& ...>& fPsis,
    const vectorField& fCutAreas,
    const std::tuple<const Field<Types>& ...>& fCutPsis,
    const pointField& ps,
    const std::tuple<const Field<Types>& ...>& pPsis,
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
            return
                Tuple2<scalar, std::tuple<Types ...>>
                (
                    cVolume,
                    tupleOp(cPsis, OpScaled<scalar>(cVolume))
                );
        }
        else
        {
            return
                Tuple2<scalar, std::tuple<Types ...>>
                (
                    scalar(0),
                    std::make_tuple(pTraits<Types>::zero ...)
                );
        }
    }

    // Averages
    const point cPAvg = OpIndirectAverage(c)(fCentres);
    auto cPsiAvgs = tupleOp(fPsis, OpIndirectAverage(c));

    // Face contributions. We use the un-cut face's centroid as the base of the
    // pyramid formed by the cut face (see !!! below). This is potentially less
    // exact than using the cut face's centroid, but it is consistent with the
    // un-cut cell volume calculation. This means if you run this function with
    // both values of "below" then the result will exactly sum to the volume of
    // the overall cell. If we used the cut-face's centroid this would not be
    // the case.
    auto result =
        cellVolumeIntegral
        (
            c,
            cAddr,
            cPAvg,
            cPsiAvgs,
            fCutAreas,
            fCentres, // !!!
            fCutPsis
        );

    // Create readably named references to the parts of the result
    scalar& cCutsVolume = result.first();
    std::tuple<Types ...>& cCutsVolumePsis = result.second();

    // Cut contributions
    const cutPoly::CellCutValues<point> cCutPValues
    (
        c,
        cAddr,
        cCuts,
        fs,
        ps,
        pAlphas,
        isoAlpha
    );
    const point cCutPAvg = OpIterableAverage()(cCutPValues);
    auto cCutPsiValues =
        tupleOp
        (
            pPsis,
            OpCellCutValues
            (
                c,
                cAddr,
                cCuts,
                fs,
                pAlphas,
                isoAlpha
            )
        );
    auto cCutPsiAvgs = tupleOp(cCutPsiValues, OpIterableAverage());

    // This method is more exact, as it uses the true centroid of the cell
    // cut. However, to obtain that centroid we have to divide by the area
    // magnitude, so this can't be generalised to types that do not support
    // division.
    const auto cCutAreaPPsis =
        faceAreaAverage
        (
            cCutPValues,
            cCutPAvg,
            std::tuple_cat(std::make_tuple(cCutPValues), cCutPsiValues),
            std::tuple_cat(std::make_tuple(cCutPAvg), cCutPsiAvgs)
        );
    const vector& cCutArea = cCutAreaPPsis.first();
    const vector& cCutCentre = std::get<0>(cCutAreaPPsis.second());
    const auto cCutPsis = tupleTail(cCutAreaPPsis.second());

    /*
    // This method is more approximate, as it uses point averages rather
    // than centroids. This does not involve division, though, so this
    // could be used with types like polynomials that only support addition
    // and multiplication.
    const vector cCutArea =
        faceAreaIntegral
        (
            cCutPValues,
            cCutPAvg,
            std::make_tuple(),
            std::make_tuple()
        ).first();
    const point& cCutCentre = cCutPAvg;
    const auto& cCutPsis = cCutPsiAvgs;
    */

    const scalar pyrVolume =
        (below ? +1 : -1)*(cCutArea & (cCutCentre - cPAvg))/3;

    cCutsVolume += pyrVolume;
    cCutsVolumePsis =
        tupleBinaryOp
        (
            cCutsVolumePsis,
            tupleBinaryOp
            (
                tupleOp
                (
                    cCutPsis,
                    OpScaled<scalar>(scalar(3)/scalar(4)*pyrVolume)
                ),
                tupleOp
                (
                    cPsiAvgs,
                    OpScaled<scalar>(scalar(1)/scalar(4)*pyrVolume)
                ),
                BinaryOpAdd()
            ),
            BinaryOpAdd()
        );

    return result;
}


// ************************************************************************* //
