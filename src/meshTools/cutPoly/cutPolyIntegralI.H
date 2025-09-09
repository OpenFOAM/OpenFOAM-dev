/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2022-2025 OpenFOAM Foundation
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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

struct OpIndex
{
    const label i_;

    inline OpIndex(const label i)
    :
        i_(i)
    {}

    template<class Type>
    inline const Type& operator()(const List<Type>& xs) const
    {
        return xs[i_];
    }
};


struct OpBegin
{
    template<class Type>
    inline auto operator()(const Type& x) const
    {
        return x.begin();
    }
};


struct OpDereference
{
    template<class Type>
    inline auto operator()(const Type& x) const
    {
        return *x;
    }
};


struct OpNext
{
    template<class Type>
    inline auto operator()(const Type& x) const
    {
        return x.next();
    }
};


template<class ScaleType>
struct OpScaled
{
    const ScaleType s_;

    inline OpScaled(const ScaleType& s)
    :
        s_(s)
    {}

    template<class Type>
    inline auto operator()(const Type& x) const
    {
        return s_*x;
    }
};


struct OpPreInner
{
    const vector& v_;

    inline OpPreInner(const vector& v)
    :
        v_(v)
    {}

    template<class Type>
    inline auto operator()(const Type& x) const
    {
        return v_ & x;
    }
};


struct OpIndirectAverage
{
    const labelUList& is_;

    inline OpIndirectAverage(const labelUList& is)
    :
        is_(is)
    {}

    template<class Container>
    inline auto operator()(const Container& xs) const
    {
        typename Container::value_type nResult =
            pTraits<typename Container::value_type>::zero;

        forAll(is_, i)
        {
            nResult += xs[is_[i]];
        }

        return nResult/is_.size();
    }
};


struct OpIterableAverage
{
    template<class Container>
    inline auto operator()(const Container& xs) const
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
};


struct OpFaceCutValues
{
    const face& f_;
    const List<labelPair>& fCuts_;
    const scalarField& pAlphas_;
    const scalar isoAlpha_;
    const bool below_;

    inline OpFaceCutValues
    (
        const face& f,
        const List<labelPair>& fCuts,
        const scalarField& pAlphas,
        const scalar isoAlpha,
        const bool below
    )
    :
        f_(f),
        fCuts_(fCuts),
        pAlphas_(pAlphas),
        isoAlpha_(isoAlpha),
        below_(below)
    {}

    template<class Type>
    inline auto operator()(const Field<Type>& pPsis) const
    {
        return
            FaceCutValues<Type>
            (
                f_,
                fCuts_,
                pPsis,
                pAlphas_,
                isoAlpha_,
                below_
            );
    }
};


struct OpCellCutValues
{
    const cell& c_;
    const cellEdgeAddressing& cAddr_;
    const labelListList& cCuts_;
    const faceList& fs_;
    const scalarField& pAlphas_;
    const scalar isoAlpha_;

    inline OpCellCutValues
    (
        const cell& c,
        const cellEdgeAddressing& cAddr,
        const labelListList& cCuts,
        const faceList& fs,
        const scalarField& pAlphas,
        const scalar isoAlpha
    )
    :
        c_(c),
        cAddr_(cAddr),
        cCuts_(cCuts),
        fs_(fs),
        pAlphas_(pAlphas),
        isoAlpha_(isoAlpha)
    {}

    template<class Type>
    inline auto operator()(const Field<Type>& pPsis) const
    {
        return
            CellCutValues<Type>
            (
                c_,
                cAddr_,
                cCuts_,
                fs_,
                pPsis,
                pAlphas_,
                isoAlpha_
            );
    }
};


struct InPlaceOpAdvance
{
    template<class Type>
    inline void operator()(Type& x) const
    {
        ++ x;
    }
};


struct BinaryOpAdd
{
    template<class Type>
    inline auto operator()(const Type& a, const Type& b) const
    {
        return a + b;
    }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace cutPoly
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

inline Foam::vector Foam::cutPoly::faceArea
(
    const face& f,
    const point& fPAvg,
    const pointField& ps
)
{
    return
        faceAreaIntegral
        (
            FaceValues<point>(f, ps),
            fPAvg,
            std::make_tuple(),
            std::make_tuple()
        ).first();
}


inline Foam::vector Foam::cutPoly::faceArea
(
    const face& f,
    const pointField& ps
)
{
    return
        faceArea
        (
            f,
            OpIndirectAverage(f)(ps),
            ps
        );
}


template<class Type>
inline Foam::Tuple2<Foam::vector, Foam::cutPoly::AreaIntegralType<Type>>
Foam::cutPoly::faceAreaIntegral
(
    const face& f,
    const point& fPAvg,
    const Type& fPsiAvg,
    const pointField& ps,
    const Field<Type>& pPsis
)
{
    auto result =
        faceAreaIntegral
        (
            FaceValues<point>(f, ps),
            fPAvg,
            std::make_tuple(FaceValues<Type>(f, pPsis)),
            std::make_tuple(fPsiAvg)
        );

    return
        Tuple2<vector, AreaIntegralType<Type>>
        (
            result.first(),
            std::get<0>(result.second())
        );
}


template<class Type>
inline Foam::Tuple2<Foam::vector, Foam::cutPoly::AreaIntegralType<Type>>
Foam::cutPoly::faceAreaIntegral
(
    const face& f,
    const pointField& ps,
    const Field<Type>& pPsis
)
{
    return
        faceAreaIntegral
        (
            f,
            OpIndirectAverage(f)(ps),
            OpIndirectAverage(f)(pPsis),
            ps,
            pPsis
        );
}


template<class Type>
inline Foam::Tuple2<Foam::vector, Type> Foam::cutPoly::faceAreaAverage
(
    const face& f,
    const point& fPAvg,
    const Type& fPsiAvg,
    const pointField& ps,
    const Field<Type>& pPsis
)
{
    auto result =
        faceAreaAverage
        (
            FaceValues<point>(f, ps),
            fPAvg,
            std::make_tuple(FaceValues<Type>(f, pPsis)),
            std::make_tuple(fPsiAvg)
        );

    return
        Tuple2<vector, Type>
        (
            result.first(),
            std::get<0>(result.second())
        );
}


template<class Type>
inline Foam::Tuple2<Foam::vector, Type> Foam::cutPoly::faceAreaAverage
(
    const face& f,
    const pointField& ps,
    const Field<Type>& pPsis
)
{
    return
        faceAreaAverage
        (
            f,
            OpIndirectAverage(f)(ps),
            OpIndirectAverage(f)(pPsis),
            ps,
            pPsis
        );
}


inline Foam::vector Foam::cutPoly::faceCutArea
(
    const face& f,
    const vector& fArea,
    const List<labelPair>& fCuts,
    const pointField& ps,
    const scalarField& pAlphas,
    const scalar isoAlpha,
    const bool below
)
{
    return
        faceCutAreaIntegral
        (
            f,
            fArea,
            std::make_tuple(),
            fCuts,
            ps,
            std::make_tuple(),
            pAlphas,
            isoAlpha,
            below
        ).first();
}


template<class Type>
inline Foam::Tuple2<Foam::vector, Foam::cutPoly::AreaIntegralType<Type>>
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
    auto result =
        faceCutAreaIntegral<Type>
        (
            f,
            fArea,
            std::make_tuple(fPsi),
            fCuts,
            ps,
            std::forward_as_tuple(pPsis),
            pAlphas,
            isoAlpha,
            below
        );

    return
        Tuple2<vector, AreaIntegralType<Type>>
        (
            result.first(),
            std::get<0>(result.second())
        );
}


inline Foam::scalar Foam::cutPoly::cellVolume
(
    const cell& c,
    const cellEdgeAddressing& cAddr,
    const point& cPAvg,
    const vectorField& fAreas,
    const pointField& fCentres
)
{
    return
        cellVolumeIntegral
        (
            c,
            cAddr,
            cPAvg,
            std::make_tuple(),
            fAreas,
            fCentres,
            std::make_tuple()
        ).first();
}


inline Foam::scalar Foam::cutPoly::cellVolume
(
    const cell& c,
    const cellEdgeAddressing& cAddr,
    const vectorField& fAreas,
    const pointField& fCentres
)
{
    return
        cellVolume
        (
            c,
            cAddr,
            OpIndirectAverage(c)(fCentres),
            fAreas,
            fCentres
        );
}


template<class Type>
inline Foam::Tuple2<Foam::scalar, Type> Foam::cutPoly::cellVolumeIntegral
(
    const cell& c,
    const cellEdgeAddressing& cAddr,
    const point& cPAvg,
    const Type& cPsiAvg,
    const vectorField& fAreas,
    const pointField& fCentres,
    const Field<Type>& fPsis
)
{
    auto result =
        cellVolumeIntegral<Type>
        (
            c,
            cAddr,
            cPAvg,
            std::make_tuple(cPsiAvg),
            fAreas,
            fCentres,
            std::forward_as_tuple(fPsis)
        );

    return
        Tuple2<scalar, Type>
        (
            result.first(),
            std::get<0>(result.second())
        );
}


template<class Type>
inline Foam::Tuple2<Foam::scalar, Type> Foam::cutPoly::cellVolumeIntegral
(
    const cell& c,
    const cellEdgeAddressing& cAddr,
    const vectorField& fAreas,
    const pointField& fCentres,
    const Field<Type>& fPsis
)
{
    return
        cellVolumeIntegral
        (
            c,
            cAddr,
            OpIndirectAverage(c)(fCentres),
            OpIndirectAverage(c)(fPsis),
            fAreas,
            fCentres,
            fPsis
        );
}


Foam::scalar Foam::cutPoly::cellCutVolume
(
    const cell& c,
    const cellEdgeAddressing& cAddr,
    const scalar cVolume,
    const labelListList& cCuts,
    const faceUList& fs,
    const vectorField& fAreas,
    const pointField& fCentres,
    const vectorField& fCutAreas,
    const pointField& ps,
    const scalarField& pAlphas,
    const scalar isoAlpha,
    const bool below
)
{
    return
        cellCutVolumeIntegral
        (
            c,
            cAddr,
            cVolume,
            std::make_tuple(),
            cCuts,
            fs,
            fAreas,
            fCentres,
            std::make_tuple(),
            fCutAreas,
            std::make_tuple(),
            ps,
            std::make_tuple(),
            pAlphas,
            isoAlpha,
            below
        ).first();
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
    const pointField& fCentres,
    const Field<Type>& fPsis,
    const vectorField& fCutAreas,
    const Field<Type>& fCutPsis,
    const pointField& ps,
    const Field<Type>& pPsis,
    const scalarField& pAlphas,
    const scalar isoAlpha,
    const bool below
)
{
    auto result =
        cellCutVolumeIntegral
        (
            c,
            cAddr,
            cVolume,
            std::make_tuple(cPsi),
            cCuts,
            fs,
            fAreas,
            fCentres,
            std::forward_as_tuple(fPsis),
            fCutAreas,
            std::forward_as_tuple(fCutPsis),
            ps,
            std::forward_as_tuple(pPsis),
            pAlphas,
            isoAlpha,
            below
        );

    return Tuple2<scalar, Type>(result.first(), std::get<0>(result.second()));
}


// ************************************************************************* //
