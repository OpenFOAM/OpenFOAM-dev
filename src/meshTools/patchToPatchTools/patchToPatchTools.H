/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2022-2023 OpenFOAM Foundation
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
    Foam::patchToPatchTools

\*---------------------------------------------------------------------------*/

#ifndef patchToPatchTools_H
#define patchToPatchTools_H

#include "distributionMap.H"
#include "remote.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace patchToPatchTools
{

//- Determine whether this intersection is confined to a single processor or
//  not, purely from the numbers of elements in the geometry
label singleProcess(const label srcSize, const label tgtSize);

//- Turn a list of send-to-processor indices into a distribution map
autoPtr<distributionMap> constructDistributionMap
(
    const labelListList& procSendIndices
);

//- Construct local addressing from the given distribution map. The result is a
//  list of remote element indices that correspond to the elements that were
//  sent by the distribution map.
List<remote> distributeAddressing(const distributionMap& map);

//- Given a local intersection addressing, determine what elements need sending
//  to each processor
labelListList procSendIndices
(
    const labelListList& tgtLocalSrcFaces,
    const List<remote>& localTgtProcFaces
);

//- Dynamic list overload of the above
labelListList procSendIndices
(
    const List<DynamicList<label>>& tgtLocalSrcFaces,
    const List<remote>& localTgtProcFaces
);

//- Trim unused elements from a distribution map. Return the map and maps
//  between the old and new construct list.
void trimDistributionMap
(
    const boolList& oldIsUsed,
    distributionMap& map,
    labelList& oldToNew,
    labelList& newToOld
);

//- Transfer list-list b into list-list a
template<class SubListA, class SubListB>
static inline void transferListList
(
    List<SubListA>& a,
    List<SubListB>& b
);

//- Reverse distribute a list-list given the map
template<class Type>
static inline void rDistributeListList
(
    const label size,
    const distributionMap& map,
    List<List<Type>>& data
);

//- Dynamic list overload of the above
template<class Type>
static inline void rDistributeListList
(
    const label size,
    const distributionMap& map,
    List<DynamicList<Type>>& data
);

//- Map local to remote indices
List<List<remote>> localToRemote
(
    const labelListList& indices,
    const List<remote>& indexToProcIndex = NullObjectRef<List<remote>>()
);

//- Dynamic list overload of the above
List<List<remote>> localToRemote
(
    const List<DynamicList<label>>& indices,
    const List<remote>& indexToProcIndex = NullObjectRef<List<remote>>()
);

//- Reverse distribute a set of target addressing
void rDistributeTgtAddressing
(
    const label tgtSize,
    const distributionMap& tgtMap,
    const List<remote>& localSrcProcFaces,
    labelListList& tgtLocalSrcFaces
);

//- Dynamic list overload of the above
void rDistributeTgtAddressing
(
    const label tgtSize,
    const distributionMap& tgtMap,
    const List<remote>& localSrcProcFaces,
    List<DynamicList<label>>& tgtLocalSrcFaces
);

//- Interpolate with normalisation
template<class Type, class LabelList, class ScalarList>
static tmp<Field<Type>> interpolate
(
    const List<LabelList>& localOtherFaces,
    const List<ScalarList>& weights,
    const autoPtr<distributionMap>& otherMapPtr,
    const Field<Type>& otherFld
);

//- Interpolate without normalisation
template<class Type, class LabelList, class ScalarList>
static tmp<Field<Type>> interpolate
(
    const List<LabelList>& localOtherFaces,
    const List<ScalarList>& weights,
    const autoPtr<distributionMap>& otherMapPtr,
    const Field<Type>& otherFld,
    const Field<Type>& leftOverFld
);

} // End namespace patchToPatchTools
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "patchToPatchToolsTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
