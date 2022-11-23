/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2022 OpenFOAM Foundation
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

Class
    Foam::cutPolyIntegral

Description
    Functions for computing integrals over cut faces and cells

SourceFiles
    cutPolyIntegralTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef cutPolyIntegral_H
#define cutPolyIntegral_H

#include "cutPolyValue.H"

#include <tuple>
#include <type_traits>
#include <utility>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace cutPoly
{

//- Compute the face-area and face-area-integrals of the given properties over
//  a face
template<template<class> class FaceValues, class ... Types>
std::tuple<vector, typename outerProduct<vector, Types>::type ...>
faceAreaIntegral
(
    const FaceValues<point>& fPs,
    const point& fPAvg,
    const std::tuple<FaceValues<Types> ...>& fPsis,
    const std::tuple<Types ...>& fPsiAvg
);

//- Compute the face-area and face-area-integral of the given property over
//  a cut-face. ToDo: Make variadic as faceAreaIntegral.
template<class Type>
Tuple2<vector, typename outerProduct<vector, Type>::type> faceCutAreaIntegral
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

//- Compute the cell-volume and cell-volume-integral of the given property over
//  a cut-cell. ToDo: Make variadic as faceAreaIntegral.
template<class Type>
Tuple2<scalar, Type> cellCutVolumeIntegral
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
);

} // End namespace cutPoly
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "cutPolyIntegralTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
