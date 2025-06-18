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

Namespace
    Foam::cutPoly

Description
    Functions for computing integrals over cut faces and cells

SourceFiles
    cutPolyIntegralTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef cutPolyIntegral_H
#define cutPolyIntegral_H

#include "cutPolyValue.H"
#include <tuple>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace cutPoly
{

//- Definition for the type of a face-area-integral
template<class Type>
using AreaIntegralType =
    typename Foam::outerProduct<Foam::vector, Type>::type;

//- Compute the face-area and face-area-integrals of the given properties
template<template<class> class FaceValues, class ... Types>
Tuple2<vector, std::tuple<AreaIntegralType<Types> ...>> faceAreaIntegral
(
    const FaceValues<point>& fPs,
    const point& fPAvg,
    const std::tuple<FaceValues<Types> ...>& fPsis,
    const std::tuple<Types ...>& fPsiAvg
);

//- Compute the face-area
template<template<class> class FaceValues>
inline vector faceArea
(
    const FaceValues<point>& fPs,
    const point& fPAvg
);

//- Compute the face-area and face-area-integral of a given property
template<template<class> class FaceValues, class Type>
inline Tuple2<vector, AreaIntegralType<Type>> faceAreaIntegral
(
    const FaceValues<point>& fPs,
    const point& fPAvg,
    const FaceValues<Type>& fPsis,
    const Type& fPsiAvg
);

//- Compute the face-cut-area and face-cut-area-integral of the given properties
template<class ... Types>
Tuple2<vector, std::tuple<AreaIntegralType<Types> ...>> faceCutAreaIntegral
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
);

//- Compute the face-cut-area
inline vector faceCutArea
(
    const face& f,
    const vector& fArea,
    const List<labelPair>& fCuts,
    const pointField& ps,
    const scalarField& pAlphas,
    const scalar isoAlpha,
    const bool below
);

//- Compute the face-cut-area and face-cut-area-integral of the given property
template<class Type>
inline Tuple2<vector, AreaIntegralType<Type>> faceCutAreaIntegral
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
);

//- Compute the cell-volume and cell-volume-integrals of the given properties
template<class ... Types>
Tuple2<scalar, std::tuple<Types ...>> cellVolumeIntegral
(
    const cell& c,
    const cellEdgeAddressing& cAddr,
    const point& cPAvg,
    const std::tuple<Types ...>& cPsiAvgs,
    const vectorField& fAreas,
    const vectorField& fCentres,
    const std::tuple<const Field<Types>& ...>& fPsis
);

//- Compute the cell-volume
inline scalar cellVolume
(
    const cell& c,
    const cellEdgeAddressing& cAddr,
    const point& cPAvg,
    const vectorField& fAreas,
    const vectorField& fCentres
);

//- Compute the cell-volume and cell-volume-integrals of the given property
template<class Type>
inline Tuple2<scalar, Type> cellVolumeIntegral
(
    const cell& c,
    const cellEdgeAddressing& cAddr,
    const point& cPAvg,
    const Type& cPsiAvg,
    const vectorField& fAreas,
    const vectorField& fCentres,
    const Field<Type>& fPsis
);

//- Compute the cell-cut-volume and cell-cut-volume-integral
//  of the given properties
template<class ... Types>
Tuple2<scalar, std::tuple<Types ...>> cellCutVolumeIntegral
(
    const cell& c,
    const cellEdgeAddressing& cAddr,
    const scalar cVolume,
    const std::tuple<Types ...>& cPsis,
    const labelListList& cCuts,
    const faceUList& fs,
    const vectorField& fAreas,
    const vectorField& fCentres,
    const std::tuple<const Field<Types>& ...>& fPsis,
    const vectorField& fCutAreas,
    const pointField& ps,
    const std::tuple<const Field<Types>& ...>& pPsis,
    const scalarField& pAlphas,
    const scalar isoAlpha,
    const bool below
);

//- Compute the cell-cut-volume
inline scalar cellCutVolume
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
);

//- Compute the cell-cut-volume and cell-cut-volume-integral
//  of the given property
template<class Type>
inline Tuple2<scalar, Type> cellCutVolumeIntegral
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
);

} // End namespace cutPoly
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "cutPolyIntegralI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "cutPolyIntegralTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
