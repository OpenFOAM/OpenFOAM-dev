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

Description
    Merge points. See below.

SourceFiles
    mergePoints.C

\*---------------------------------------------------------------------------*/

#ifndef mergePoints_H
#define mergePoints_H

#include "scalar.H"
#include "labelList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Function mergePoints Declaration
\*---------------------------------------------------------------------------*/

//- Sorts and merges points. All points closer than/equal mergeTol get merged.
//  Returns the number of unique points and a map from old to new.
template<class Type>
label mergePoints
(
    const UList<Type>& points,
    const scalar mergeTol,
    const bool verbose,
    labelList& pointMap,
    const Type& origin = Type::zero
);

//- Sorts and merges points. Determines new points. Returns true if anything
//  merged (though newPoints still sorted even if not merged).
template<class Type>
bool mergePoints
(
    const UList<Type>& points,
    const scalar mergeTol,
    const bool verbose,
    labelList& pointMap,
    List<Type>& newPoints,
    const Type& origin = Type::zero
);

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "mergePoints.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
