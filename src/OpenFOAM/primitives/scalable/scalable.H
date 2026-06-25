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

Struct
    Foam::scalable

Description
    Trait to identify types which are "scalable"; i.e., that can be
    multiply-equals-d with a scalar. This has to be done explicitly, rather
    than by identifying the existence of the operator with decltype, in order
    to exclude bool and int and similar.

\*---------------------------------------------------------------------------*/

#ifndef scalable_H
#define scalable_H

#include <type_traits>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class scalable Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
struct scalable : public std::false_type {};

template<class Type, class Void = void>
using enableIfScalable = std::enable_if_t<scalable<Type>::value, Void>;

template<class Type, class Void = void>
using enableIfNotScalable = std::enable_if_t<!scalable<Type>::value, Void>;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
