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

template<template<class> class FaceValues>
inline Foam::vector Foam::cutPoly::faceArea
(
    const FaceValues<point>& fPs,
    const point& fPAvg
)
{
    return
        faceAreaIntegral
        (
            fPs,
            fPAvg,
            std::make_tuple(),
            std::make_tuple()
        ).first();
}


template<template<class> class FaceValues, class Type>
inline Foam::Tuple2<Foam::vector, Foam::cutPoly::AreaIntegralType<Type>>
Foam::cutPoly::faceAreaIntegral
(
    const FaceValues<point>& fPs,
    const point& fPAvg,
    const FaceValues<Type>& fPsis,
    const Type& fPsiAvg
)
{
    auto result =
        faceAreaIntegral
        (
            fPs,
            fPAvg,
            std::make_tuple(fPsis),
            std::make_tuple(fPsiAvg)
        );

    return
        Tuple2<vector, AreaIntegralType<Type>>
        (
            result.first(),
            std::get<0>(result.second())
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
    const vectorField& fCentres
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


template<class Type>
inline Foam::Tuple2<Foam::scalar, Type> Foam::cutPoly::cellVolumeIntegral
(
    const cell& c,
    const cellEdgeAddressing& cAddr,
    const point& cPAvg,
    const Type& cPsiAvg,
    const vectorField& fAreas,
    const vectorField& fCentres,
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


Foam::scalar Foam::cutPoly::cellCutVolume
(
    const cell& c,
    const cellEdgeAddressing& cAddr,
    const scalar cVolume,
    const labelListList& cCuts,
    const faceUList& fs,
    const vectorField& fAreas,
    const vectorField& fCentres,
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
    const vectorField& fCentres,
    const vectorField& fPsis,
    const vectorField& fCutAreas,
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
            ps,
            std::forward_as_tuple(pPsis),
            pAlphas,
            isoAlpha,
            below
        );

    return Tuple2<scalar, Type>(result.first(), std::get<0>(result.second()));
}


// ************************************************************************* //
