/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2025 OpenFOAM Foundation
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

#include "IOmanip.H"
#include "indexedOctree.H"
#include "labelIOField.H"
#include "LagrangianMesh.H"
#include "LagrangianMeshLocation.H"
#include "LagrangianModels.H"
#include "ListOps.H"
#include "meshObjects.H"
#include "Time.H"
#include "tracking.H"
#include "treeDataCell.H"
#include "debug.H"

#include "internalLagrangianPatch.H"
#include "nonConformalCyclicLagrangianPatch.H"
#include "nonConformalProcessorCyclicLagrangianPatch.H"
#include "nonConformalErrorLagrangianPatch.H"

#include "calculatedLagrangianPatchFields.H"
#include "internalLagrangianFieldSources.H"
#include "zeroLagrangianFieldSources.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(LagrangianMesh, 0);

    const word LagrangianMesh::prefix("Lagrangian");

    const word LagrangianMesh::coordinatesName("coordinates");
    const word LagrangianMesh::positionName("position");
    const word LagrangianMesh::stateName("state");
    const word LagrangianMesh::fractionName("fraction");

    template<>
    const char* NamedEnum<LagrangianMesh::permutationAlgorithm, 2>::names[] =
        {"copy", "inPlace"};

    const NamedEnum<LagrangianMesh::permutationAlgorithm, 2>
        LagrangianMesh::permutationAlgorithmNames_;

    LagrangianMesh::permutationAlgorithm
        LagrangianMesh::permutationAlgorithm_ =
        Foam::debug::namedEnumOptimisationSwitch
        (
            (LagrangianMesh::typeName + "Permutation").c_str(),
            LagrangianMesh::permutationAlgorithmNames_,
            LagrangianMesh::permutationAlgorithm::copy
        );

    template<>
    const char* NamedEnum<LagrangianMesh::partitioningAlgorithm, 2>::names[] =
        {"bin", "quick"};

    const NamedEnum<LagrangianMesh::partitioningAlgorithm, 2>
        LagrangianMesh::partitioningAlgorithmNames_;

    LagrangianMesh::partitioningAlgorithm
        LagrangianMesh::partitioningAlgorithm_ =
        Foam::debug::namedEnumOptimisationSwitch
        (
            (LagrangianMesh::typeName + "Partitioning").c_str(),
            LagrangianMesh::partitioningAlgorithmNames_,
            LagrangianMesh::partitioningAlgorithm::bin
        );
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::LagrangianMesh::printGroups(const bool header) const
{
    checkPtr(offsetsPtr_, "Offsets");
    const labelList& offsets = offsetsPtr_();

    // Shorthand casts of group indices
    static const label completei =
        static_cast<label>(LagrangianGroup::complete);
    static const label inInternalMeshi =
        static_cast<label>(LagrangianGroup::inInternalMesh);
    static const label onPatchZeroi =
        static_cast<label>(LagrangianGroup::onPatchZero);

    // Determine the number of non-processor patches
    label nNonProcPatches = 0;
    for (; nNonProcPatches < boundary().size(); ++ nNonProcPatches)
    {
        const LagrangianPatch& p = boundary()[nNonProcPatches];

        if (isA<processorLagrangianPatch>(p)) break;
    }

    // Determine if there are processor patches and/or if any patches are
    // referred to by processor-cyclics
    bool hasProcPatches = false;
    boolList patchIsReferred(nNonProcPatches, false);
    for (label patchi = nNonProcPatches; patchi < boundary().size(); ++ patchi)
    {
        const LagrangianPatch& p = boundary()[patchi];

        if (isA<processorCyclicLagrangianPatch>(p))
        {
            patchIsReferred
            [
                refCast<const processorCyclicLagrangianPatch>(p)
               .referPatchIndex()
            ] = true;
        }
        else // isA<processorLagrangianPatch>(p)
        {
            hasProcPatches = true;
        }
    }
    reduce(hasProcPatches, orOp<bool>());
    Pstream::listCombineGather(patchIsReferred, orEqOp<bool>());
    Pstream::listCombineScatter(patchIsReferred);

    // Get the indices of patches that are referred to by processor-cyclics
    const labelList referredPatches(findIndices(patchIsReferred, true));
    const labelList patchReferred(invert(nNonProcPatches, referredPatches));

    // Determine the indices of the groups in the table after those relating to
    // the global patches
    const label onProcPatchi = onPatchZeroi + nNonProcPatches;
    const label onProcCyclicPatchZeroi = onProcPatchi + hasProcPatches;
    const label toBeRemovedi = onProcCyclicPatchZeroi + referredPatches.size();

    // Build the column names
    wordList columnNames(toBeRemovedi + 1);
    columnNames[completei] = Foam::name(LagrangianGroup::complete);
    columnNames[inInternalMeshi] = Foam::name(LagrangianGroup::inInternalMesh);
    for (label patchi = 0; patchi < nNonProcPatches; ++ patchi)
    {
        columnNames[onPatchZeroi + patchi] =
            isA<nonConformalErrorLagrangianPatch>(boundary()[patchi])
         || isA<internalLagrangianPatch>(boundary()[patchi])
          ? word::null
          : boundary()[patchi].name();
    }
    if (hasProcPatches)
    {
        columnNames[onProcPatchi] = "(processor)";
    }
    forAll(referredPatches, referredPatchi)
    {
        columnNames[onProcCyclicPatchZeroi + referredPatchi] =
            "(processor)"
          + boundary()[referredPatches[referredPatchi]].name();
    }
    columnNames[toBeRemovedi] = Foam::name(LagrangianGroup::toBeRemoved);

    // Print the column names if the header is requested, and finish
    if (header)
    {
        forAll(columnNames, i)
        {
            if (i && columnNames[i].size()) Info<< ' ';

            Info<< columnNames[i];
        }

        return;
    }

    // Construct the numbers of elements in each column
    labelList columnNumbers(toBeRemovedi + 1, 0);
    columnNumbers[completei] = offsets[completei + 1] - offsets[completei];
    columnNumbers[inInternalMeshi] =
        offsets[inInternalMeshi + 1] - offsets[inInternalMeshi];
    for (label patchi = 0; patchi < nNonProcPatches; ++ patchi)
    {
        columnNumbers[onPatchZeroi + patchi] =
            offsets[onPatchZeroi + patchi + 1] - offsets[onPatchZeroi + patchi];
    }
    for (label patchi = nNonProcPatches; patchi < boundary().size(); ++ patchi)
    {
        const LagrangianPatch& p = boundary()[patchi];

        if (isA<processorCyclicLagrangianPatch>(p))
        {
            columnNumbers
            [
                patchReferred
                [
                    refCast<const processorCyclicLagrangianPatch>(p)
                   .referPatchIndex()
                ]
            ] +=
                offsets[onPatchZeroi + patchi + 1]
              - offsets[onPatchZeroi + patchi];
        }
        else // isA<processorLagrangianPatch>(p)
        {
            columnNumbers[onProcPatchi] +=
                offsets[onPatchZeroi + patchi + 1]
              - offsets[onPatchZeroi + patchi];
        }
    }
    Pstream::listCombineGather(columnNumbers, plusEqOp<label>());
    Pstream::listCombineScatter(columnNumbers);

    // Print the numbers
    forAll(columnNames, i)
    {
        if (i && columnNames[i].size()) Info<< ' ';

        const string::size_type columnNumberDigits =
            Foam::name(columnNumbers[i]).size();

        if (columnNumberDigits > columnNames[i].size())
        {
            Info << string(columnNames[i].size(), '#').c_str();
        }
        else
        {
            Info << setw(columnNames[i].size()) << columnNumbers[i];
        }
    }
}


Foam::labelList Foam::LagrangianMesh::partitionBin
(
    labelList& offsets,
    const List<LagrangianState>& states
) const
{
    // Set the zero index so we only partition from the first previously
    // incomplete element onwards
    const label i0 = offsets[1];

    // Sum the numbers of elements in each group and store in the offsets array
    offsets = 0;
    for (label i = i0; i < states.size(); ++ i)
    {
        const label groupi = stateToGroupi(states[i]);
        offsets[groupi + 1] ++;
    }

    // Cumulative sum, starting at i0, to generate the offsets for the
    // non-complete subsets of elements
    offsets[0] = i0;
    for (label groupi = 0; groupi < nGroups(); ++ groupi)
    {
        offsets[groupi + 1] += offsets[groupi];
    }

    // Insert each element into the permutation. Increment the offsets to keep
    // track of the current insertion position within each group.
    labelList permutation(states.size() - i0);
    for (label i = i0; i < states.size(); ++ i)
    {
        const label groupi = stateToGroupi(states[i]);
        permutation[offsets[groupi] - i0] = i;
        ++ offsets[groupi];
    }

    // The offsets are now shifted by one. Copy them back a place and set the
    // first offset to zero, rather than i0, so that the offsets refer to the
    // entire mesh, not just the non-complete subsets.
    for (label groupi = nGroups(); groupi > 0; -- groupi)
    {
        offsets[groupi] = offsets[groupi - 1];
    }
    offsets[0] = 0;

    // Return the permutation
    return permutation;
}


Foam::labelList Foam::LagrangianMesh::partitionQuick
(
    labelList& offsets,
    const List<LagrangianState>& states
) const
{
    // Set the zero index so we only partition from the first previously
    // incomplete element onwards
    const label i0 = offsets[1];

    // Initialise the offsets, excluding all previously complete elements
    offsets = -1;
    offsets[0] = i0;
    offsets[nGroups()] = states.size();

    // Initialise a permutation
    labelList permutation(states.size() - i0);
    forAll(permutation, i)
    {
        permutation[i] = i + i0;
    }

    // Compute the number of iterations needed to order the groups
    label nIterations = 0;
    for (label groupi = 1; groupi < nGroups(); groupi *= 2)
    {
        ++ nIterations;
    }

    // Partition
    label nDivisions = 1;
    for (label iterationi = 0; iterationi < nIterations; ++ iterationi)
    {
        for (label divisioni = 0; divisioni < nDivisions; ++ divisioni)
        {
            const label pivot = (2*divisioni + 1)*nGroups()/(2*nDivisions);
            if (offsets[pivot] != -1)
            {
                continue;
            }

            const label pivot0 = divisioni*nGroups()/nDivisions;
            const label pivot1 = (divisioni + 1)*nGroups()/nDivisions;
            label i = offsets[pivot0] - i0, j = offsets[pivot1] - i0;
            while (i < j)
            {
                while
                (
                    i < j
                 && stateToGroupi(states[permutation[i]]) < pivot
                )
                {
                    ++ i;
                }
                while
                (
                    i < j
                 && stateToGroupi(states[permutation[j - 1]]) >= pivot
                )
                {
                    -- j;
                }
                if (i < j)
                {
                    Swap(permutation[i], permutation[j - 1]);
                }
            }

            offsets[pivot] = i + i0;
        }

        nDivisions *= 2;
    }

    // Re-include previously complete elements into the complete group
    offsets[0] = 0;

    // Return the permutation
    return permutation;
}


void Foam::LagrangianMesh::permuteAndResizeFields(const labelList& permutation)
{
    wordHashSet permutedFieldNames;
    #define PERMUTE_TYPE_FIELDS(Type, GeoField)                                \
    {                                                                          \
        HashTable<GeoField<Type>*> fields                                      \
        (                                                                      \
            lookupCurrentFields<GeoField<Type>>()                              \
        );                                                                     \
                                                                               \
        forAllIter(typename HashTable<GeoField<Type>*>, fields, iter)          \
        {                                                                      \
            if (permutedFieldNames.found(iter()->name())) continue;            \
                                                                               \
            permutedFieldNames.insert(iter()->name());                         \
                                                                               \
            permuteList(permutation, iter()->primitiveFieldRef());             \
                                                                               \
            resizeContainer(iter()->primitiveFieldRef());                      \
        }                                                                      \
    }
    PERMUTE_TYPE_FIELDS(label, LagrangianField);
    FOR_ALL_FIELD_TYPES(PERMUTE_TYPE_FIELDS, LagrangianField);
    PERMUTE_TYPE_FIELDS(label, LagrangianDynamicField);
    FOR_ALL_FIELD_TYPES(PERMUTE_TYPE_FIELDS, LagrangianDynamicField);
    PERMUTE_TYPE_FIELDS(label, LagrangianInternalField);
    FOR_ALL_FIELD_TYPES(PERMUTE_TYPE_FIELDS, LagrangianInternalField);
    PERMUTE_TYPE_FIELDS(label, LagrangianInternalDynamicField);
    FOR_ALL_FIELD_TYPES(PERMUTE_TYPE_FIELDS, LagrangianInternalDynamicField);
    #undef PERMUTE_TYPE_FIELDS
}


template<class Type>
void Foam::LagrangianMesh::permuteList
(
    const labelList& permutation,
    UList<Type>& list
)
{
    switch (permutationAlgorithm_)
    {
        case permutationAlgorithm::copy:
            permuteListCopy(permutation, list);
            break;

        case permutationAlgorithm::inPlace:
            permuteListInPlace(permutation, list);
            break;
    }
}


template<class Type>
void Foam::LagrangianMesh::permuteListCopy
(
    const labelList& permutation,
    UList<Type>& list
)
{
    if (permutation.empty()) return;

    const label i0 = list.size() - permutation.size();

    SubList<Type>(list, permutation.size(), i0) =
        List<Type>(UIndirectList<Type>(list, permutation));
}


template<class Type>
void Foam::LagrangianMesh::permuteListInPlace
(
    const labelList& permutation,
    UList<Type>& list
)
{
    if (permutation.empty()) return;

    const label i0 = list.size() - permutation.size();

    labelList& permutationRef = const_cast<labelList&>(permutation);

    label i = 0;
    Type t = list[i + i0];

    forAll(permutation, permutationi)
    {
        bool end = permutationRef[i] < 0;

        while (permutationRef[i] < 0) ++ i;

        if (end) t = list[i + i0];

        //// Forward
        //Swap(list[permutationRef[i]], t);

        // Reverse
        list[i + i0] =
            permutationRef[permutationRef[i] - i0] < 0
          ? t : list[permutationRef[i]];

        permutationRef[i] = - permutationRef[i] - 1;
        i = - permutationRef[i] - 1 - i0;
    }

    forAll(permutationRef, permutationi)
    {
        permutationRef[permutationi] = - permutationRef[permutationi] - 1;
    }
}


template<class Container>
void Foam::LagrangianMesh::resizeContainer(Container& container) const
{
    container.resize(offsetsPtr_->last());
}


Foam::LagrangianSubMesh Foam::LagrangianMesh::append
(
    const barycentricField& coordinates,
    const labelField& celli,
    const labelField& facei,
    const labelField& faceTrii
)
{
    clearPosition();

    const LagrangianSubMesh appendMesh
    (
        *this,
        LagrangianGroup::none,
        coordinates.size(),
        size()
    );

    if (statesPtr_.valid())
    {
        states().resize(appendMesh.end(), LagrangianState::none);
    }

    if (receivePatchFacePtr_.valid())
    {
        receivePatchFacePtr_().resize(appendMesh.end(), label(-1));
    }

    if (receivePositionPtr_.valid())
    {
        receivePositionPtr_().resize(appendMesh.end(), point::nan);
    }

    coordinates_.append(coordinates);
    celli_.append(celli);
    facei_.append(facei);
    faceTrii_.append(faceTrii);

    subAll_.size_ = size();

    return appendMesh;
}


Foam::LagrangianSubMesh Foam::LagrangianMesh::append
(
    const labelList& parents
)
{
    clearPosition();

    const LagrangianSubMesh appendMesh
    (
        *this,
        LagrangianGroup::none,
        parents.size(),
        size()
    );

    if (statesPtr_.valid())
    {
        states().resize(appendMesh.end());
        appendMesh.sub(static_cast<List<LagrangianState>&>(states())) =
            UIndirectList<LagrangianState>(states(), parents)();
    }

    if (receivePatchFacePtr_.valid())
    {
        receivePatchFacePtr_().resize(appendMesh.end());
        appendMesh.sub(static_cast<List<label>&>(receivePatchFacePtr_())) =
            UIndirectList<label>(receivePatchFacePtr_(), parents)();
    }

    if (receivePositionPtr_.valid())
    {
        receivePositionPtr_().resize(appendMesh.end());
        appendMesh.sub(static_cast<List<point>&>(receivePositionPtr_())) =
            UIndirectList<point>(receivePositionPtr_(), parents)();
    }

    coordinates_.resize(appendMesh.end());
    appendMesh.sub(static_cast<List<barycentric>&>(coordinates_)) =
        UIndirectList<barycentric>(coordinates_, parents)();
    celli_.resize(appendMesh.end());
    appendMesh.sub(static_cast<labelList&>(celli_)) =
        UIndirectList<label>(celli_, parents)();
    facei_.resize(appendMesh.end());
    appendMesh.sub(static_cast<labelList&>(facei_)) =
        UIndirectList<label>(facei_, parents)();
    faceTrii_.resize(appendMesh.end());
    appendMesh.sub(static_cast<labelList&>(faceTrii_)) =
        UIndirectList<label>(faceTrii_, parents)();

    subAll_.size_ = size();

    return appendMesh;
}


Foam::wordHashSet Foam::LagrangianMesh::appendSpecifiedFields
(
    const LagrangianSubMesh&
) const
{
    return wordHashSet();
}


void Foam::LagrangianMesh::injectUnspecifiedFields
(
    const LagrangianInjection& injection,
    const LagrangianSubMesh& injectionMesh,
    const wordHashSet& specifiedFieldNames
)
{
    // Inject values for unspecified fields using their source conditions
    wordHashSet injectedFieldNames;
    #define INJECT_TYPE_FIELDS(Type, GeoField)                                 \
    {                                                                          \
        HashTable<GeoField<Type>*> fields                                      \
        (                                                                      \
            lookupCurrentFields<GeoField<Type>>()                              \
        );                                                                     \
                                                                               \
        forAllIter(typename HashTable<GeoField<Type>*>, fields, iter)          \
        {                                                                      \
            if (specifiedFieldNames.found(iter()->name())) continue;           \
                                                                               \
            injectedFieldNames.insert(iter()->name());                         \
                                                                               \
            /* Resize the field */                                             \
            iter()->resize(injectionMesh.end());                               \
                                                                               \
            /* Use the source condition to set the new values */               \
            injectionMesh.sub(*iter()).ref() =                                 \
                iter()->sources()[injection.name()].value                      \
                (                                                              \
                    injection,                                                 \
                    injectionMesh                                              \
                );                                                             \
        }                                                                      \
    }
    INJECT_TYPE_FIELDS(label, LagrangianField);
    FOR_ALL_FIELD_TYPES(INJECT_TYPE_FIELDS, LagrangianField);
    INJECT_TYPE_FIELDS(label, LagrangianDynamicField);
    FOR_ALL_FIELD_TYPES(INJECT_TYPE_FIELDS, LagrangianDynamicField);
    #undef INJECT_TYPE_FIELDS

    // Special handling for a retained state field
    #define INJECT_STATE_FIELD(GeoField)                                       \
    {                                                                          \
        if (foundObject<GeoField<label>>(stateName))                           \
        {                                                                      \
            GeoField<label>& state =                                           \
                lookupObjectRef<GeoField<label>>(stateName);                   \
                                                                               \
            injectedFieldNames.insert(stateName);                              \
                                                                               \
            /* Resize the field */                                             \
            state.resize(injectionMesh.end());                                 \
                                                                               \
            /* Set unknown state */                                            \
            injectionMesh.sub(state).ref() =                                   \
                static_cast<label>(LagrangianState::none);                     \
        }                                                                      \
    }
    INJECT_STATE_FIELD(LagrangianInternalField);
    INJECT_STATE_FIELD(LagrangianInternalDynamicField);
    #undef INJECT_STATE_FIELD

    // Make a table of internal fields for which values could not be set
    wordHashSet internalFieldNames;
    #define INSERT_INTERNAL_FIELD_NAMES(Type, GeoField)                        \
    {                                                                          \
        HashTable<GeoField<Type>*> fields                                      \
        (                                                                      \
            lookupCurrentFields<GeoField<Type>>()                              \
        );                                                                     \
                                                                               \
        forAllIter(typename HashTable<GeoField<Type>*>, fields, iter)          \
        {                                                                      \
            if (specifiedFieldNames.found(iter()->name())) continue;           \
            if (injectedFieldNames.found(iter()->name())) continue;            \
                                                                               \
            internalFieldNames.insert(iter()->name());                         \
        }                                                                      \
    }
    INSERT_INTERNAL_FIELD_NAMES(label, LagrangianInternalField);
    FOR_ALL_FIELD_TYPES(INSERT_INTERNAL_FIELD_NAMES, LagrangianInternalField);
    INSERT_INTERNAL_FIELD_NAMES(label, LagrangianInternalDynamicField);
    FOR_ALL_FIELD_TYPES
    (
        INSERT_INTERNAL_FIELD_NAMES,
        LagrangianInternalDynamicField
    );
    #undef INSERT_INTERNAL_FIELD_NAMES

    if (!internalFieldNames.empty())
    {
        FatalErrorInFunction
            << "Internal fields " << internalFieldNames.sortedToc() << " are "
            << "registered during an injection event. These fields do not "
            << "contain source conditions so their new values cannot be "
            << "assigned. Either specify these fields' values in the calling "
            << "code, or ensure that they are not registered."
            << exit(FatalError);
    }
}


void Foam::LagrangianMesh::injectUnspecifiedFields
(
    const LagrangianSubMesh& injectionMesh,
    const wordHashSet& specifiedFieldNames
)
{
    wordHashSet fieldNames;
    #define INSERT_FIELD_NAMES(Type, GeoField) \
        fieldNames.insert(lookupCurrentFields<GeoField<Type>>().toc());
    INSERT_FIELD_NAMES(label, LagrangianField);
    FOR_ALL_FIELD_TYPES(INSERT_FIELD_NAMES, LagrangianField);
    INSERT_FIELD_NAMES(label, LagrangianDynamicField);
    FOR_ALL_FIELD_TYPES(INSERT_FIELD_NAMES, LagrangianDynamicField);
    INSERT_FIELD_NAMES(label, LagrangianInternalField);
    FOR_ALL_FIELD_TYPES(INSERT_FIELD_NAMES, LagrangianInternalField);
    INSERT_FIELD_NAMES(label, LagrangianInternalDynamicField);
    FOR_ALL_FIELD_TYPES(INSERT_FIELD_NAMES, LagrangianInternalDynamicField);
    #undef INSERT_FIELD_NAMES

    if (!fieldNames.empty())
    {
        FatalErrorInFunction
            << "Fields " << fieldNames.sortedToc() << " are registered during "
            << "an injection event. These fields' new values cannot be "
            << "assigned. Either specify these fields' values in the calling "
            << "code, or ensure that they are not registered."
            << exit(FatalError);
    }
}


void Foam::LagrangianMesh::birthUnspecifiedFields
(
    const labelList& parents,
    const LagrangianSubMesh& birthMesh,
    const wordHashSet& specifiedFieldNames
)
{
    wordHashSet birthedFieldNames;
    #define BIRTH_TYPE_FIELDS(Type, GeoField)                                  \
    {                                                                          \
        HashTable<GeoField<Type>*> fields                                      \
        (                                                                      \
            lookupCurrentFields<GeoField<Type>>()                              \
        );                                                                     \
                                                                               \
        forAllIter(typename HashTable<GeoField<Type>*>, fields, iter)          \
        {                                                                      \
            if (specifiedFieldNames.found(iter()->name())) continue;           \
            if (birthedFieldNames.found(iter()->name())) continue;             \
                                                                               \
            birthedFieldNames.insert(iter()->name());                          \
                                                                               \
            /* Resize the field */                                             \
            iter()->resize(birthMesh.end());                                   \
                                                                               \
            /* Map values from the parent elements */                          \
            birthMesh.sub(*iter()).ref().primitiveFieldRef() =                 \
                Field<Type>(UIndirectList<Type>(*iter(), parents)());          \
        }                                                                      \
    }
    BIRTH_TYPE_FIELDS(label, LagrangianField);
    FOR_ALL_FIELD_TYPES(BIRTH_TYPE_FIELDS, LagrangianField);
    BIRTH_TYPE_FIELDS(label, LagrangianDynamicField);
    FOR_ALL_FIELD_TYPES(BIRTH_TYPE_FIELDS, LagrangianDynamicField);
    BIRTH_TYPE_FIELDS(label, LagrangianInternalField);
    FOR_ALL_FIELD_TYPES(BIRTH_TYPE_FIELDS, LagrangianInternalField);
    BIRTH_TYPE_FIELDS(label, LagrangianInternalDynamicField);
    FOR_ALL_FIELD_TYPES(BIRTH_TYPE_FIELDS, LagrangianInternalDynamicField);
    #undef BIRTH_TYPE_FIELDS
}


void Foam::LagrangianMesh::changer::constructNonConformal() const
{
    bool haveNccPatches = false;

    forAll(mesh_.boundary(), patchi)
    {
        const polyPatch& pp = mesh_.mesh().boundaryMesh()[patchi];

        if (isA<nonConformalCyclicPolyPatch>(pp))
        {
            haveNccPatches = true;
            break;
        }
    }

    if (!haveNccPatches) return;

    mesh_.origPatchNccPatchisPtr_.set
    (
        new labelListList
        (
            mesh_.boundary().size(),
            labelList()
        )
    );
    mesh_.origPatchNccPatchesPtr_.set
    (
        new List<UPtrList<const nonConformalCyclicPolyPatch>>
        (
            mesh_.boundary().size(),
            UPtrList<const nonConformalCyclicPolyPatch>()
        )
    );
    mesh_.nccPatchProcNccPatchisPtr_.set
    (
        new labelListList
        (
            mesh_.boundary().size(),
            labelList(Pstream::nProcs(), -1)
        )
    );

    forAll(mesh_.boundary(), patchi)
    {
        const polyPatch& pp = mesh_.mesh().boundaryMesh()[patchi];

        if (isA<nonConformalCyclicPolyPatch>(pp))
        {
            const nonConformalCyclicPolyPatch& nccPp =
                refCast<const nonConformalCyclicPolyPatch>(pp);

            mesh_.origPatchNccPatchisPtr_()
                [nccPp.origPatchIndex()].append(patchi);
            mesh_.origPatchNccPatchesPtr_()
                [nccPp.origPatchIndex()].append(&nccPp);

            mesh_.nccPatchProcNccPatchisPtr_()
                [nccPp.index()][Pstream::myProcNo()] =
                nccPp.index();

            if (nccPp.owner()) nccPp.rays();
        }
        else if (isA<nonConformalProcessorCyclicPolyPatch>(pp))
        {
            const nonConformalProcessorCyclicPolyPatch& ncpcPp =
                refCast<const nonConformalProcessorCyclicPolyPatch>(pp);

            mesh_.nccPatchProcNccPatchisPtr_()
                [ncpcPp.referPatchIndex()][ncpcPp.neighbProcNo()] =
                ncpcPp.index();
        }
    }

    mesh_.receivePatchFacePtr_.set
    (
        new DynamicList<label>(mesh_.size(), label(-1))
    );

    mesh_.receivePositionPtr_.set
    (
        new DynamicList<point>(mesh_.size(), point::nan)
    );
}


void Foam::LagrangianMesh::changer::constructBehind() const
{
    mesh_.fractionBehindPtr_.set
    (
        new LagrangianDynamicField<scalar>
        (
            IOobject
            (
                fractionName + "Behind",
                mesh_.time().name(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensioned<scalar>(dimless, scalar(0)),
            wordList
            (
                mesh_.boundary().size(),
                calculatedLagrangianPatchScalarField::typeName
            ),
            wordList::null(),
            LagrangianModels::New(mesh_).modelTypeFieldSourceTypes
            <
                LagrangianInjection,
                zeroLagrangianScalarFieldSource,
                LagrangianSource,
                internalLagrangianScalarFieldSource
            >()
        )
    );

    mesh_.nTracksBehindPtr_.set
    (
        new LagrangianDynamicField<label>
        (
            IOobject
            (
                "nTracksBehind",
                mesh_.time().name(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensioned<label>(dimless, label(0)),
            wordList
            (
                mesh_.boundary().size(),
                calculatedLagrangianPatchLabelField::typeName
            ),
            wordList::null(),
            LagrangianModels::New(mesh_).modelTypeFieldSourceTypes
            <
                LagrangianInjection,
                zeroLagrangianLabelFieldSource,
                LagrangianSource,
                internalLagrangianLabelFieldSource
            >()
        )
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::LagrangianMesh::LagrangianMesh
(
    const polyMesh& mesh,
    const word& name,
    const IOobject::readOption readOption,
    const IOobject::writeOption writeOption
)
:
    LagrangianMesh
    (
        mesh,
        name,
        mesh.boundaryMesh().types(),
        readOption,
        writeOption
    )
{}


Foam::LagrangianMesh::LagrangianMesh
(
    const polyMesh& mesh,
    const word& name,
    const wordList& wantedPatchTypes,
    const IOobject::readOption readOption,
    const IOobject::writeOption writeOption
)
:
    objectRegistry
    (
        IOobject
        (
            name,
            mesh.time().name(),
            prefix,
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        )
    ),
    GeoMesh<polyMesh>(mesh),
    mesh_(mesh),
    coordinates_
    (
        IOobject
        (
            coordinatesName,
            time().name(),
            *this,
            readOption,
            IOobject::AUTO_WRITE
        )
    ),
    celli_
    (
        IOobject
        (
            "cell",
            time().name(),
            *this,
            readOption,
            IOobject::AUTO_WRITE
        )
    ),
    facei_
    (
        IOobject
        (
            "face",
            time().name(),
            *this,
            readOption,
            IOobject::AUTO_WRITE
        )
    ),
    faceTrii_
    (
        IOobject
        (
            "faceTri",
            time().name(),
            *this,
            readOption,
            IOobject::AUTO_WRITE
        )
    ),
    boundary_(*this, mesh.boundaryMesh(), wantedPatchTypes),
    subAll_
    (
        LagrangianSubMesh
        (
            *this,
            LagrangianGroup::none,
            size(),
            label(0),
            label(0)
        )
    ),
    statesPtr_(nullptr),
    offsetsPtr_(nullptr),
    subMeshIndex_(0),
    schemesPtr_(nullptr)
{
    writeOpt() = writeOption;

    checkFieldSize(coordinates_);
    checkFieldSize(celli_);

    // Read cellFace file if face does not exist. This is useful for testing,
    // as it is much easier to write a cellFace file by hand.
    if (!celli_.empty() && facei_.empty())
    {
        labelIOField cellFacei
        (
            IOobject
            (
                "cellFace",
                time().name(),
                *this,
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            )
        );

        if (!cellFacei.empty())
        {
            checkFieldSize(cellFacei);
            facei_.resize(cellFacei.size());
            forAll(facei_, i)
            {
                facei_[i] = mesh_.cells()[celli_[i]][cellFacei[i]];
            }
        }
    }

    checkFieldSize(facei_);
    checkFieldSize(faceTrii_);

    // Ask for the tetBasePtIs and oldCellCentres to trigger all processors to
    // build them, otherwise, if some processors have no elements then there is
    // a comms mismatch.
    mesh_.tetBasePtIs();
    mesh_.oldCellCentres();
}


Foam::LagrangianMesh::changer::changer
(
    LagrangianMesh& mesh,
    const LagrangianState state
)
:
    mesh_(mesh)
{
    mesh_.statesPtr_.set
    (
        new DynamicList<LagrangianState>(mesh_.size(), state)
    );

    mesh_.offsetsPtr_.set
    (
        new labelList(mesh_.nGroups() + 1, 0)
    );

    for
    (
        label groupi = mesh_.stateToGroupi(state) + 1;
        groupi < mesh_.nGroups() + 1;
        ++ groupi
    )
    {
        mesh_.offsetsPtr_()[groupi] = mesh_.size();
    }

    constructNonConformal();

    mesh_.subMeshIndex_ = 0;

    Info<< indent;
    mesh_.printGroups(true);
    Info<< endl << indent;
    mesh_.printGroups(false);
    Info<< endl;

    forAll(mesh_.boundary(), patchi)
    {
        mesh_.boundary()[patchi].partition();
    }

    constructBehind();
}


Foam::LagrangianMesh::changer::changer
(
    LagrangianMesh& mesh,
    const List<LagrangianState>& states
)
:
    mesh_(mesh)
{
    mesh_.statesPtr_.set
    (
        new DynamicList<LagrangianState>(states)
    );

    mesh_.offsetsPtr_.set
    (
        new labelList(mesh_.nGroups() + 1, 0)
    );

    constructNonConformal();

    mesh_.subMeshIndex_ = 0;

    Info<< indent;
    mesh_.printGroups(true);
    Info<< endl;
    mesh_.partition();

    constructBehind();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::LagrangianMesh::~LagrangianMesh()
{}


Foam::LagrangianMesh::changer::~changer()
{
    mesh_.statesPtr_.clear();
    mesh_.offsetsPtr_.clear();

    mesh_.origPatchNccPatchisPtr_.clear();
    mesh_.origPatchNccPatchesPtr_.clear();
    mesh_.nccPatchProcNccPatchisPtr_.clear();
    mesh_.receivePatchFacePtr_.clear();
    mesh_.receivePositionPtr_.clear();

    mesh_.fractionBehindPtr_.clear();
    mesh_.nTracksBehindPtr_.clear();

    mesh_.subMeshIndex_ = 0;

    forAll(mesh_.boundary(), patchi)
    {
        mesh_.boundary()[patchi].partition();
    }
}


Foam::LagrangianMesh::linearDisplacement::linearDisplacement
(
    const LagrangianSubVectorField& linear
)
:
    linear_(linear)
{}


Foam::LagrangianMesh::parabolicDisplacement::parabolicDisplacement
(
    const LagrangianSubVectorField& linear,
    const LagrangianSubVectorField& quadratic
)
:
    linear_(linear),
    quadratic_(quadratic)
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

const Foam::LagrangianSchemes& Foam::LagrangianMesh::schemes() const
{
    if (!schemesPtr_.valid())
    {
        schemesPtr_ = new LagrangianSchemes(*this);
    }

    return schemesPtr_();
}


const Foam::LagrangianSolution& Foam::LagrangianMesh::solution() const
{
    if (!solutionPtr_.valid())
    {
        solutionPtr_ = new LagrangianSolution(*this);
    }

    return solutionPtr_();
}


Foam::labelList Foam::LagrangianMesh::subMeshGlobalSizes() const
{
    // Determine the number of global groups. Assumes processor patches come
    // after all global patches.
    label nGlobalGroups = static_cast<label>(LagrangianGroup::onPatchZero);
    forAll(boundary(), patchi)
    {
        const polyPatch& pp = boundary()[patchi].patch();
        if (isA<processorPolyPatch>(pp)) break;
        nGlobalGroups ++;
    }

    // Extend by one for the to-be-removed group
    nGlobalGroups ++;

    static const label completei =
        static_cast<label>(LagrangianGroup::complete);
    static const label inInternalMeshi =
        static_cast<label>(LagrangianGroup::inInternalMesh);
    static const label onPatchZeroi =
        static_cast<label>(LagrangianGroup::onPatchZero);

    // Create a list of sizes of the global groups
    labelList sizes(nGlobalGroups, 0);
    sizes[completei] = sub(LagrangianGroup::complete).size();
    sizes[inInternalMeshi] = sub(LagrangianGroup::inInternalMesh).size();
    forAll(boundary(), patchi)
    {
        const polyPatch& pp = boundary()[patchi].patch();
        if (isA<processorPolyPatch>(pp)) break;
        sizes[onPatchZeroi + patchi] = boundary()[patchi].mesh().size();
    }
    // (to-be-removed group last)
    sizes.last() =
        size() - boundary()[nGlobalGroups - 2 - onPatchZeroi].mesh().start();

    // Sum over all processes
    Pstream::listCombineGather(sizes, plusEqOp<label>());
    Pstream::listCombineScatter(sizes);

    // Construct the result. This includes processor patch groups, on which a
    // "size" of -1 is set.
    labelList result(nGroups(), -1);
    for (label i = 0; i < nGlobalGroups - 1; ++ i)
    {
        result[i] = sizes[i];
    }
    // (to-be-removed group last)
    result.last() = sizes.last();

    return result;
}


Foam::tmp<Foam::LagrangianVectorInternalField>
Foam::LagrangianMesh::position() const
{
    tmp<LagrangianVectorInternalField> tresult =
        LagrangianVectorInternalField::New
        (
            positionName,
            *this,
            dimLength
        );
    LagrangianVectorInternalField& result = tresult.ref();

    forAll(coordinates_, i)
    {
        result[i] = position(i);
    }

    return tresult;
}


Foam::point Foam::LagrangianMesh::position(const label i) const
{
    return
        tracking::position
        (
            mesh_,
            coordinates_[i], celli_[i], facei_[i], faceTrii_[i], 1
        );
}


Foam::LagrangianMesh::location Foam::LagrangianMesh::locate
(
    const point& position,
    barycentric& coordinates,
    label& celli,
    label& facei,
    label& faceTrii,
    const scalar fraction
) const
{
    // Look for a containing cell and set the process if found
    remote procCelli;
    procCelli.elementi = mesh().cellTree().findInside(position);
    procCelli.proci = procCelli.elementi >= 0 ? Pstream::myProcNo() : -1;

    // Pick a unique processor
    reduce(procCelli, remote::firstProcOp());

    // Find the tetrahedron and the local coordinates
    location result = location::outsideMesh;
    if (procCelli.proci == Pstream::myProcNo())
    {
        result =
            tracking::locate
            (
                mesh_, position,
                coordinates, celli, facei, faceTrii, fraction
            )
          ? location::inCell
          : location::onBoundary;
    }

    // Communicate the location flag and return
    reduce(result, LagrangianMeshLocation::closestOp());

    return result;
}


Foam::List<Foam::LagrangianMesh::location> Foam::LagrangianMesh::locate
(
    const List<point>& position,
    List<barycentric>& coordinates,
    labelList& celli,
    labelList& facei,
    labelList& faceTrii,
    const scalarList& fraction
) const
{
    // Look for containing cells and set the process if found
    List<remote> procCelli(position.size());
    forAll(position, i)
    {
        procCelli[i].elementi = mesh().cellTree().findInside(position[0]);
        procCelli[i].proci =
            procCelli[i].elementi >= 0 ? Pstream::myProcNo() : -1;
    }

    // Pick unique processors
    reduce(procCelli, ListOp<remote::firstProcOp>());

    // Find the tetrahedra and the local coordinates
    List<location> result(position.size(), location::outsideMesh);
    forAll(position, i)
    {
        if (procCelli[i].proci != Pstream::myProcNo()) continue;

        result =
            tracking::locate
            (
                mesh_, position[i],
                coordinates[i], celli[i], facei[i], faceTrii[i], fraction[i]
            )
          ? location::inCell
          : location::onBoundary;
    }

    // Communicate the location flags and return
    reduce(result, ListOp<LagrangianMeshLocation::closestOp>());

    return result;
}


void Foam::LagrangianMesh::partition()
{
    clearPosition();

    // Get the offsets
    checkPtr(offsetsPtr_, "Offsets");
    labelList& offsets = offsetsPtr_();

    // Partition the state offsets and states and create a permutation
    labelList permutation;
    switch (partitioningAlgorithm_)
    {
        case partitioningAlgorithm::bin:
            permutation = partitionBin(offsets, states());
            break;

        case partitioningAlgorithm::quick:
            permutation = partitionQuick(offsets, states());
            break;
    }

    // Print the updated states
    Info<< indent;
    printGroups(false);
    Info<< endl;

    // Apply the permutation to the states and positions
    permuteList(permutation, states());
    permuteList(permutation, coordinates_);
    permuteList(permutation, celli_);
    permuteList(permutation, facei_);
    permuteList(permutation, faceTrii_);

    // Apply the permutation to the non-conformal receive information (if any)
    if (receivePatchFacePtr_.valid())
    {
        permuteList(permutation, receivePatchFacePtr_());
    }
    if (receivePositionPtr_.valid())
    {
        permuteList(permutation, receivePositionPtr_());
    }

    // Check that the states are ordered
    #ifdef FULLDEBUG
    for (label groupi = 0; groupi < nGroups(); ++ groupi)
    {
        for (label i = offsets[groupi]; i < offsets[groupi + 1]; ++ i)
        {
            if (stateToGroupi(states()[i]) != groupi)
            {
                FatalErrorInFunction
                    << "Partitioning failed"
                    << exit(FatalError);
            }
        }
    }
    #endif

    // Reclaim space by removing the toBeRemoved group
    offsets[nGroups()] = offsets[nGroups() - 1];

    // Resize the states and positions
    resizeContainer(states());
    resizeContainer(coordinates_);
    resizeContainer(celli_);
    resizeContainer(facei_);
    resizeContainer(faceTrii_);

    // Resize the non-conformal receive information (if any)
    if (receivePatchFacePtr_.valid())
    {
        resizeContainer(receivePatchFacePtr_());
    }
    if (receivePositionPtr_.valid())
    {
        resizeContainer(receivePositionPtr_());
    }

    // Permute and resize the fields
    permuteAndResizeFields(permutation);

    // Update the patches
    forAll(boundary(), patchi)
    {
        boundary()[patchi].partition();
    }

    // Update the sub-all mesh
    subAll_.size_ = size();
}


template<class Displacement>
void Foam::LagrangianMesh::track
(
    const List<LagrangianState>& endState,
    const Displacement& displacement,
    const LagrangianSubScalarField& deltaFraction,
    LagrangianSubScalarSubField& fraction
)
{
    clearPosition();

    // The fraction is about to change. Ensure the previous values are stored
    // to facilitate subsequent calculations.
    fraction.oldTime();

    // Track each element in the sub-mesh in turn
    forAll(fraction, subi)
    {
        const label i = subi + fraction.mesh().start();

        // Track to completion or the next face
        Tuple2<bool, scalar> onFaceAndF =
            tracking::toFace
            (
                mesh_, displacement(subi), deltaFraction[subi],
                coordinates_[i], celli_[i], facei_[i], faceTrii_[i],
                fraction[subi],
                fractionBehindPtr_()[i], nTracksBehindPtr_()[i],
                debug
              ? static_cast<const string&>(name() + " #" + Foam::name(i))
              : NullObjectRef<string>()
            );

        // Update the state
        if (!onFaceAndF.first())
        {
            states()[i] = endState[subi];
        }
        else if (mesh_.isInternalFace(facei_[i]))
        {
            states()[i] = LagrangianState::onInternalFace;
        }
        else // if (<on a boundary face>)
        {
            // Determine the index of the patch that was tracked to
            label patchi =
                mesh_.boundaryMesh().patchIndices()
                [
                    facei_[i] - mesh_.nInternalFaces()
                ];

            // If this patch has non-conformal cyclics associated with it, then
            // search through them and see if any was hit. If we find one that
            // does, override the patch index variable.
            if
            (
                origPatchNccPatchisPtr_.valid()
             && origPatchNccPatchisPtr_()[patchi].size()
            )
            {
                // Get the current position
                const point sendPosition =
                    tracking::position
                    (
                        mesh_,
                        coordinates_[i], celli_[i], facei_[i], faceTrii_[i],
                        fraction[subi]
                    );

                // Get the displacement of the location that was hit
                const vector sendDisplacement =
                    tracking::faceNormalAndDisplacement
                    (
                        mesh_,
                        coordinates_[i], celli_[i], facei_[i], faceTrii_[i],
                        fraction[subi]
                    ).second();

                // Use ray searching on each non-conformal cyclic in turn
                forAll(origPatchNccPatchisPtr_()[patchi], patchNccPatchi)
                {
                    const label nccPatchi =
                        origPatchNccPatchisPtr_()[patchi][patchNccPatchi];
                    const nonConformalCyclicPolyPatch& nccPp =
                        origPatchNccPatchesPtr_()[patchi][patchNccPatchi];

                    point receivePosition;
                    const remote receiveProcAndFace =
                        nccPp.ray
                        (
                            fraction[subi],
                            nccPp.origPatch().whichFace(facei_[i]),
                            sendPosition,
                            displacement(subi, onFaceAndF.second())
                          - fraction[subi]*sendDisplacement,
                            receivePosition
                        );

                    const label receiveProci = receiveProcAndFace.proci;

                    if (receiveProci == -1) continue;

                    const label receiveFacei = receiveProcAndFace.elementi;

                    receivePatchFacePtr_()[i] = receiveFacei;
                    receivePositionPtr_()[i] = receivePosition;

                    patchi =
                        nccPatchProcNccPatchisPtr_()[nccPatchi][receiveProci];

                    break;
                }
            }

            // Set the state to that of the identified patch
            states()[i] =
                static_cast<LagrangianState>
                (
                    static_cast<label>(LagrangianState::onPatchZero)
                  + patchi
                );
        }
    }
}


template
void Foam::LagrangianMesh::track<Foam::LagrangianMesh::linearDisplacement>
(
    const List<LagrangianState>& endState,
    const linearDisplacement& displacement,
    const LagrangianSubScalarField& deltaFraction,
    LagrangianSubScalarSubField& fraction
);


template
void Foam::LagrangianMesh::track<Foam::LagrangianMesh::parabolicDisplacement>
(
    const List<LagrangianState>& endState,
    const parabolicDisplacement& displacement,
    const LagrangianSubScalarField& deltaFraction,
    LagrangianSubScalarSubField& fraction
);


void Foam::LagrangianMesh::crossFaces
(
    const LagrangianScalarInternalDynamicField& fraction
)
{
    clearPosition();

    // Internal face-crossings
    const LagrangianSubMesh incompleteMesh
    (
        sub(LagrangianGroup::inInternalMesh)
    );

    forAll(incompleteMesh, subi)
    {
        const label i = subi + incompleteMesh.start();

        if (states()[i] != LagrangianState::onInternalFace) continue;

        // Cross the face
        tracking::crossInternalFace
        (
            mesh_,
            coordinates_[i], celli_[i], facei_[i], faceTrii_[i]
        );

        // Update the state
        states()[i] = LagrangianState::inCell;
    }

    // Patch-face crossings and boundary condition evaluations
    if
    (
        Pstream::defaultCommsType == Pstream::commsTypes::blocking
     || Pstream::defaultCommsType == Pstream::commsTypes::nonBlocking
    )
    {
        PstreamBuffers pBufs(Pstream::defaultCommsType);

        // Initialise patch crossings
        forAll(boundary(), patchi)
        {
            boundary()[patchi].initEvaluate(pBufs, *this, fraction);
        }

        // Initialise patch field evaluations
        #define INIT_EVAL_TYPE_PATCH_FIELDS(Type, GeoField)                    \
        {                                                                      \
            HashTable<GeoField<Type>*> fields                                  \
            (                                                                  \
                lookupCurrentFields<GeoField<Type>>()                          \
            );                                                                 \
                                                                               \
            forAllIter(typename HashTable<GeoField<Type>*>, fields, iter)      \
            {                                                                  \
                forAll(boundary(), patchi)                                     \
                {                                                              \
                    iter()->boundaryFieldRef()[patchi].initEvaluate            \
                    (                                                          \
                        pBufs,                                                 \
                        fraction                                               \
                    );                                                         \
                }                                                              \
            }                                                                  \
        }
        INIT_EVAL_TYPE_PATCH_FIELDS(label, LagrangianField);
        FOR_ALL_FIELD_TYPES(INIT_EVAL_TYPE_PATCH_FIELDS, LagrangianField);
        INIT_EVAL_TYPE_PATCH_FIELDS(label, LagrangianDynamicField);
        FOR_ALL_FIELD_TYPES(INIT_EVAL_TYPE_PATCH_FIELDS,LagrangianDynamicField);
        #undef INIT_EVAL_TYPE_PATCH_FIELDS

        // Block for any outstanding requests
        pBufs.finishedSends();

        // Patch crossings
        forAll(boundary(), patchi)
        {
            boundary()[patchi].evaluate(pBufs, *this, fraction);
        }

        // Patch field evaluations
        #define EVAL_TYPE_PATCH_FIELDS(Type, GeoField)                         \
        {                                                                      \
            HashTable<GeoField<Type>*> fields                                  \
            (                                                                  \
                lookupCurrentFields<GeoField<Type>>()                          \
            );                                                                 \
                                                                               \
            forAllIter(typename HashTable<GeoField<Type>*>, fields, iter)      \
            {                                                                  \
                forAll(boundary(), patchi)                                     \
                {                                                              \
                    iter()->boundaryFieldRef()[patchi].evaluate                \
                    (                                                          \
                        pBufs,                                                 \
                        fraction                                               \
                    );                                                         \
                }                                                              \
            }                                                                  \
        }
        EVAL_TYPE_PATCH_FIELDS(label, LagrangianField);
        FOR_ALL_FIELD_TYPES(EVAL_TYPE_PATCH_FIELDS, LagrangianField);
        EVAL_TYPE_PATCH_FIELDS(label, LagrangianDynamicField);
        FOR_ALL_FIELD_TYPES(EVAL_TYPE_PATCH_FIELDS, LagrangianDynamicField);
        #undef EVAL_TYPE_PATCH_FIELDS
    }
    else
    {
        FatalErrorInFunction
            << "Unsupported communications type "
            << Pstream::commsTypeNames[Pstream::defaultCommsType]
            << exit(FatalError);
    }
}


void Foam::LagrangianMesh::reset(const bool initial, const bool final)
{
    clearPosition();

    if (initial && final) return;

    if (initial)
    {
        coordinates_.storeOldTimes();
        coordinates_.oldTime();
        celli_.storeOldTimes();
        celli_.oldTime();
        facei_.storeOldTimes();
        facei_.oldTime();
        faceTrii_.storeOldTimes();
        faceTrii_.oldTime();

        #define OLD_TIME_TYPE_FIELDS(Type, GeoField)                           \
        {                                                                      \
            HashTable<GeoField<Type>*> fields                                  \
            (                                                                  \
                lookupCurrentFields<GeoField<Type>>()                          \
            );                                                                 \
                                                                               \
            forAllIter(typename HashTable<GeoField<Type>*>, fields, iter)      \
            {                                                                  \
                iter()->storeOldTimes();                                       \
                iter()->oldTime();                                             \
            }                                                                  \
        }
        OLD_TIME_TYPE_FIELDS(label, LagrangianDynamicField);
        FOR_ALL_FIELD_TYPES(OLD_TIME_TYPE_FIELDS, LagrangianDynamicField);
        #undef OLD_TIME_TYPE_FIELDS
    }

    if (!initial)
    {
        static_cast<DynamicField<barycentric>&>(coordinates_) =
            coordinates_.oldTime();
        static_cast<DynamicField<label>&>(celli_) = celli_.oldTime();
        static_cast<DynamicField<label>&>(facei_) = facei_.oldTime();
        static_cast<DynamicField<label>&>(faceTrii_) = faceTrii_.oldTime();

        #define RESET_OLD_TIME_TYPE_FIELDS(Type, GeoField)                     \
        {                                                                      \
            HashTable<GeoField<Type>*> fields                                  \
            (                                                                  \
                lookupCurrentFields<GeoField<Type>>()                          \
            );                                                                 \
                                                                               \
            forAllIter(typename HashTable<GeoField<Type>*>, fields, iter)      \
            {                                                                  \
                iter()->reset(iter()->oldTime());                              \
            }                                                                  \
        }
        RESET_OLD_TIME_TYPE_FIELDS(label, LagrangianDynamicField);
        FOR_ALL_FIELD_TYPES(RESET_OLD_TIME_TYPE_FIELDS, LagrangianDynamicField);
        #undef RESET_OLD_TIME_TYPE_FIELDS
    }

    if (final)
    {
        coordinates_.clearOldTimes();
        celli_.clearOldTimes();
        facei_.clearOldTimes();
        faceTrii_.clearOldTimes();

        #define CLEAR_OLD_TIME_TYPE_FIELDS(Type, GeoField)                     \
        {                                                                      \
            HashTable<GeoField<Type>*> fields                                  \
            (                                                                  \
                lookupCurrentFields<GeoField<Type>>()                          \
            );                                                                 \
                                                                               \
            forAllIter(typename HashTable<GeoField<Type>*>, fields, iter)      \
            {                                                                  \
                iter()->clearOldTimes();                                       \
            }                                                                  \
        }
        CLEAR_OLD_TIME_TYPE_FIELDS(label, LagrangianDynamicField);
        FOR_ALL_FIELD_TYPES(CLEAR_OLD_TIME_TYPE_FIELDS, LagrangianDynamicField);
        #undef CLEAR_OLD_TIME_TYPE_FIELDS
    }

    wordHashSet fieldNames;
    #define INSERT_FIELD_NAMES(Type, GeoField) \
        fieldNames.insert(lookupCurrentFields<GeoField<Type>>(true).toc());
    INSERT_FIELD_NAMES(label, LagrangianField);
    FOR_ALL_FIELD_TYPES(INSERT_FIELD_NAMES, LagrangianField);
    INSERT_FIELD_NAMES(label, LagrangianInternalField);
    FOR_ALL_FIELD_TYPES(INSERT_FIELD_NAMES, LagrangianInternalField);
    INSERT_FIELD_NAMES(label, LagrangianInternalDynamicField);
    FOR_ALL_FIELD_TYPES(INSERT_FIELD_NAMES, LagrangianInternalDynamicField);
    #undef INSERT_FIELD_NAMES

    if (!fieldNames.empty())
    {
        FatalErrorInFunction
            << "Non-dynamic and/or internal fields " << fieldNames.sortedToc()
            << " are registered during a reset. Only fields relating to "
            << "fundamental state should be present at this time, and these "
            << "should all be dynamic non-internal fields" << exit(FatalError);
    }

    subAll_.size_ = size();
}


void Foam::LagrangianMesh::clear()
{
    clearPosition();

    if (statesPtr_.valid())
    {
        states().clear();
    }

    coordinates_.clear();
    celli_.clear();
    facei_.clear();
    faceTrii_.clear();

    subAll_.size_ = 0;
}


Foam::labelList Foam::LagrangianMesh::partition
(
    const label nGroups,
    const UList<labelPair>& elementsGroups
)
{
    if (elementsGroups.empty()) return labelList(nGroups + 1, size());

    // !!! Not implemented as yet. This is needed for models which
    // instantaneously modify or remove elements. The only instantaneous models
    // implemented thus far are injection models which only create elements.
    // This will need implementing for breakup models and similar.
    NotImplemented;
    return labelList();
}


void Foam::LagrangianMesh::remove(const UList<label>& elements)
{
    List<labelPair> elementsGroups(elements.size());
    forAll(elements, i)
    {
        elementsGroups[i] = labelPair(elements[i], 0);
    }

    const labelList offsets = partition(1, elementsGroups);

    remove(offsets[1] - offsets[0]);
}


void Foam::LagrangianMesh::remove(const label nElements)
{
    if (nElements == 0) return;

    // !!! Not implemented as yet. This is needed for models which
    // instantaneously remove elements. The only instantaneous models
    // implemented thus far are injection models which only create elements.
    // This will need implementing for breakup models and similar.
    NotImplemented;
}


void Foam::LagrangianMesh::clearPosition()
{
    if (foundObject<LagrangianVectorInternalField>(positionName))
    {
        lookupObjectRef<LagrangianVectorInternalField>(positionName).checkOut();
    }
}


void Foam::LagrangianMesh::storePosition()
{
    if (!this->foundObject<LagrangianVectorInternalField>(positionName))
    {
        this->position().ptr()->store();
    }
}


void Foam::LagrangianMesh::topoChange(const polyTopoChangeMap& map)
{
    NotImplemented;

    meshObjects::topoChange<LagrangianMesh>(*this, map);
}


void Foam::LagrangianMesh::mapMesh(const polyMeshMap& map)
{
    NotImplemented;

    meshObjects::mapMesh<LagrangianMesh>(*this, map);
}


void Foam::LagrangianMesh::distribute(const polyDistributionMap& map)
{
    NotImplemented;

    meshObjects::distribute<LagrangianMesh>(*this, map);
}


bool Foam::LagrangianMesh::writeObject
(
    IOstream::streamFormat fmt,
    IOstream::versionNumber ver,
    IOstream::compressionType cmp,
    const bool write
) const
{
    return objectRegistry::writeObject(fmt, ver, cmp, write);
}


bool Foam::LagrangianMesh::write(const bool write) const
{
    return objectRegistry::write(write);
}


// ************************************************************************* //
