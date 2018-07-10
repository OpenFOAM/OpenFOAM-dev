/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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

InClass
    Foam

Description
    Spatial transformation functions for primitive fields.

SourceFiles
    transformList.C

\*---------------------------------------------------------------------------*/

#ifndef transformList_H
#define transformList_H

#include "transform.H"
#include "List.H"
#include "Map.H"
#include "EdgeMap.H"
#include "tensorField.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- Extend transform to work on list.
template<class T>
List<T> transform
(
    const tensor& rotTensor,
    const UList<T>& field
);

//- Apply transformation to list. Either single transformation tensor
//  or one tensor per element.
template<class T>
void transformList(const tensor&, UList<T>&);
template<class T>
void transformList(const tensorField&, UList<T>&);

template<class T>
void transformList(const tensor&, Map<T>&);
template<class T>
void transformList(const tensorField&, Map<T>&);

template<class T>
void transformList(const tensor&, EdgeMap<T>&);
template<class T>
void transformList(const tensorField&, EdgeMap<T>&);


// Specialisations for bool
template<>
inline void transformList(const tensor&, UList<bool>&)
{}
template<>
inline void transformList(const tensorField&, UList<bool>&)
{}
template<>
inline void transformList(const tensor&, Map<bool>&)
{}
template<>
inline void transformList(const tensorField&, Map<bool>&)
{}
template<>
inline void transformList(const tensor&, EdgeMap<bool>&)
{}
template<>
inline void transformList(const tensorField&, EdgeMap<bool>&)
{}


// Specialisations for label
template<>
inline void transformList(const tensor&, labelUList&)
{}
template<>
inline void transformList(const tensorField&, labelUList&)
{}
template<>
inline void transformList(const tensor&, Map<label>&)
{}
template<>
inline void transformList(const tensorField&, Map<label>&)
{}
template<>
inline void transformList(const tensor&, EdgeMap<label>&)
{}
template<>
inline void transformList(const tensorField&, EdgeMap<label>&)
{}


// Specialisations for scalar
template<>
inline void transformList(const tensor&, UList<scalar>&)
{}
template<>
inline void transformList(const tensorField&, UList<scalar>&)
{}
template<>
inline void transformList(const tensor&, Map<scalar>&)
{}
template<>
inline void transformList(const tensorField&, Map<scalar>&)
{}
template<>
inline void transformList(const tensor&, EdgeMap<scalar>&)
{}
template<>
inline void transformList(const tensorField&, EdgeMap<scalar>&)
{}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "transformList.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
