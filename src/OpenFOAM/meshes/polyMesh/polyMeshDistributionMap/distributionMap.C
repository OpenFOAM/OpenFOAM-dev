/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2022 OpenFOAM Foundation
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

#include "distributionMap.H"
#include "globalIndexAndTransform.H"
#include "transformField.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(distributionMap, 0);
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<>
void Foam::distributionMap::transform::operator()
(
    const transformer&,
    const bool,
    List<label>&
) const
{}
template<>
void Foam::distributionMap::transform::operator()
(
    const coupledPolyPatch&,
    UList<label>&
) const
{}
template<>
void Foam::distributionMap::transform::operator()
(
    const coupledPolyPatch&,
    Map<label>&
) const
{}
template<>
void Foam::distributionMap::transform::operator()
(
    const coupledPolyPatch&,
    EdgeMap<label>&
) const
{}


template<>
void Foam::distributionMap::transform::operator()
(
    const transformer&,
    const bool,
    List<scalar>&
) const
{}
template<>
void Foam::distributionMap::transform::operator()
(
    const coupledPolyPatch&,
    UList<scalar>&
) const
{}
template<>
void Foam::distributionMap::transform::operator()
(
    const coupledPolyPatch&,
    Map<scalar>&
) const
{}
template<>
void Foam::distributionMap::transform::operator()
(
    const coupledPolyPatch&,
    EdgeMap<scalar>&
) const
{}


template<>
void Foam::distributionMap::transform::operator()
(
    const transformer&,
    const bool,
    List<bool>&
) const
{}
template<>
void Foam::distributionMap::transform::operator()
(
    const coupledPolyPatch&,
    UList<bool>&
) const
{}
template<>
void Foam::distributionMap::transform::operator()
(
    const coupledPolyPatch&,
    Map<bool>&
) const
{}
template<>
void Foam::distributionMap::transform::operator()
(
    const coupledPolyPatch&,
    EdgeMap<bool>&
) const
{}


void Foam::distributionMap::printLayout(Ostream& os) const
{
    distributionMapBase::printLayout(os);

    forAll(transformElements_, trafoI)
    {
        if (transformElements_[trafoI].size() > 0)
        {
            os  << "transform " << trafoI << ':' << endl
                << "    start : " << transformStart_[trafoI] << endl
                << "    size  : " << transformElements_[trafoI].size() << endl;
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::distributionMap::distributionMap()
:
    distributionMapBase()
{}


Foam::distributionMap::distributionMap
(
    const label constructSize,
    labelListList&& subMap,
    labelListList&& constructMap,
    const bool subHasFlip,
    const bool constructHasFlip
)
:
    distributionMapBase
    (
        constructSize,
        move(subMap),
        move(constructMap),
        subHasFlip,
        constructHasFlip
    )
{}


Foam::distributionMap::distributionMap
(
    const label constructSize,
    labelListList&& subMap,
    labelListList&& constructMap,
    labelListList&& transformElements,
    labelList&& transformStart,
    const bool subHasFlip,
    const bool constructHasFlip
)
:
    distributionMapBase
    (
        constructSize,
        move(subMap),
        move(constructMap),
        subHasFlip,
        constructHasFlip
    ),
    transformElements_(move(transformElements)),
    transformStart_(move(transformStart))
{}


Foam::distributionMap::distributionMap
(
    const labelList& sendProcs,
    const labelList& recvProcs
)
:
    distributionMapBase(sendProcs, recvProcs)
{}


Foam::distributionMap::distributionMap
(
    const globalIndex& globalNumbering,
    labelList& elements,
    List<Map<label>>& compactMap,
    const int tag
)
:
    distributionMapBase
    (
        globalNumbering,
        elements,
        compactMap,
        tag
    )
{}


Foam::distributionMap::distributionMap
(
    const globalIndex& globalNumbering,
    labelListList& cellCells,
    List<Map<label>>& compactMap,
    const int tag
)
:
    distributionMapBase
    (
        globalNumbering,
        cellCells,
        compactMap,
        tag
    )
{}


Foam::distributionMap::distributionMap
(
    const globalIndex& globalNumbering,
    labelList& elements,
    const globalIndexAndTransform& globalTransforms,
    const labelPairList& transformedElements,
    labelList& transformedIndices,
    List<Map<label>>& compactMap,
    const int tag
)
:
    distributionMapBase()
{
    // Construct per processor compact addressing of the global elements
    // needed. The ones from the local processor are not included since
    // these are always all needed.
    calcCompactAddressing
    (
        globalNumbering,
        elements,
        compactMap
    );

    // Add all (non-local) transformed elements needed.
    forAll(transformedElements, i)
    {
        labelPair elem = transformedElements[i];
        label proci = globalTransforms.processor(elem);
        if (proci != Pstream::myProcNo())
        {
            label index = globalTransforms.index(elem);
            label nCompact = compactMap[proci].size();
            compactMap[proci].insert(index, nCompact);
        }
    }


    // Exchange what I need with processor that supplies it. Renumber elements
    // into compact numbering
    labelList compactStart;
    exchangeAddressing
    (
        tag,
        globalNumbering,
        elements,
        compactMap,
        compactStart
    );


    // Renumber the transformed elements
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Count per transformIndex
    label nTrafo = globalTransforms.transformPermutations().size();
    labelList nPerTransform(nTrafo, 0);
    forAll(transformedElements, i)
    {
        labelPair elem = transformedElements[i];
        label trafoI = globalTransforms.transformIndex(elem);
        nPerTransform[trafoI]++;
    }
    // Offset per transformIndex
    transformStart_.setSize(nTrafo);
    transformElements_.setSize(nTrafo);
    forAll(transformStart_, trafoI)
    {
        transformStart_[trafoI] = constructSize_;
        constructSize_ += nPerTransform[trafoI];
        transformElements_[trafoI].setSize(nPerTransform[trafoI]);
    }

    // Sort transformed elements into their new slot.
    nPerTransform = 0;

    transformedIndices.setSize(transformedElements.size());
    forAll(transformedElements, i)
    {
        labelPair elem = transformedElements[i];
        label proci = globalTransforms.processor(elem);
        label index = globalTransforms.index(elem);
        label trafoI = globalTransforms.transformIndex(elem);

        // Get compact index for untransformed element
        label rawElemI =
        (
            proci == Pstream::myProcNo()
          ? index
          : compactMap[proci][index]
        );

        label& n = nPerTransform[trafoI];
        // index of element to transform
        transformElements_[trafoI][n] = rawElemI;
        // destination of transformed element
        transformedIndices[i] = transformStart_[trafoI]+n;
        n++;
    }

    if (debug)
    {
        printLayout(Pout);
    }
}


Foam::distributionMap::distributionMap
(
    const globalIndex& globalNumbering,
    labelListList& cellCells,
    const globalIndexAndTransform& globalTransforms,
    const List<labelPairList>& transformedElements,
    labelListList& transformedIndices,
    List<Map<label>>& compactMap,
    const int tag
)
:
    distributionMapBase()
{
    // Construct per processor compact addressing of the global elements
    // needed. The ones from the local processor are not included since
    // these are always all needed.
    calcCompactAddressing
    (
        globalNumbering,
        cellCells,
        compactMap
    );

    // Add all (non-local) transformed elements needed.
    forAll(transformedElements, celli)
    {
        const labelPairList& elems = transformedElements[celli];

        forAll(elems, i)
        {
            label proci = globalTransforms.processor(elems[i]);
            if (proci != Pstream::myProcNo())
            {
                label index = globalTransforms.index(elems[i]);
                label nCompact = compactMap[proci].size();
                compactMap[proci].insert(index, nCompact);
            }
        }
    }


    // Exchange what I need with processor that supplies it. Renumber elements
    // into compact numbering
    labelList compactStart;
    exchangeAddressing
    (
        tag,
        globalNumbering,
        cellCells,
        compactMap,
        compactStart
    );


    // Renumber the transformed elements
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Count per transformIndex
    label nTrafo = globalTransforms.transformPermutations().size();
    labelList nPerTransform(nTrafo, 0);
    forAll(transformedElements, celli)
    {
        const labelPairList& elems = transformedElements[celli];

        forAll(elems, i)
        {
            label trafoI = globalTransforms.transformIndex(elems[i]);
            nPerTransform[trafoI]++;
        }
    }
    // Offset per transformIndex
    transformStart_.setSize(nTrafo);
    transformElements_.setSize(nTrafo);
    forAll(transformStart_, trafoI)
    {
        transformStart_[trafoI] = constructSize_;
        constructSize_ += nPerTransform[trafoI];
        transformElements_[trafoI].setSize(nPerTransform[trafoI]);
    }

    // Sort transformed elements into their new slot.
    nPerTransform = 0;

    transformedIndices.setSize(transformedElements.size());
    forAll(transformedElements, celli)
    {
        const labelPairList& elems = transformedElements[celli];
        transformedIndices[celli].setSize(elems.size());

        forAll(elems, i)
        {
            label proci = globalTransforms.processor(elems[i]);
            label index = globalTransforms.index(elems[i]);
            label trafoI = globalTransforms.transformIndex(elems[i]);

            // Get compact index for untransformed element
            label rawElemI =
            (
                proci == Pstream::myProcNo()
              ? index
              : compactMap[proci][index]
            );

            label& n = nPerTransform[trafoI];
            // index of element to transform
            transformElements_[trafoI][n] = rawElemI;
            // destination of transformed element
            transformedIndices[celli][i] = transformStart_[trafoI]+n;
            n++;
        }
    }

    if (debug)
    {
        printLayout(Pout);
    }
}


Foam::distributionMap::distributionMap(const distributionMap& map)
:
    distributionMapBase(map),
    transformElements_(map.transformElements_),
    transformStart_(map.transformStart_)
{}


Foam::distributionMap::distributionMap(distributionMap&& map)
:
    distributionMapBase(move(map)),
    transformElements_(move(map.transformElements_)),
    transformStart_(move(map.transformStart_))
{}


Foam::distributionMap::distributionMap(Istream& is)
{
    is  >> *this;
}


Foam::autoPtr<Foam::distributionMap> Foam::distributionMap::clone() const
{
    return autoPtr<distributionMap>(new distributionMap(*this));
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::label Foam::distributionMap::whichTransform(const label index)
const
{
    return findLower(transformStart_, index+1);
}


void Foam::distributionMap::transfer(distributionMap& rhs)
{
    distributionMapBase::transfer(rhs);
    transformElements_.transfer(rhs.transformElements_);
    transformStart_.transfer(rhs.transformStart_);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

void Foam::distributionMap::operator=(const distributionMap& rhs)
{
    // Check for assignment to self
    if (this == &rhs)
    {
        FatalErrorInFunction
            << "Attempted assignment to self"
            << abort(FatalError);
    }
    distributionMapBase::operator=(rhs);
    transformElements_ = rhs.transformElements_;
    transformStart_ = rhs.transformStart_;
}


// * * * * * * * * * * * * * * Istream Operator  * * * * * * * * * * * * * * //

Foam::Istream& Foam::operator>>(Istream& is, distributionMap& map)
{
    is.fatalCheck("operator>>(Istream&, distributionMap&)");

    is  >> static_cast<distributionMapBase&>(map)
        >> map.transformElements_ >> map.transformStart_;

    return is;
}


// * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<(Ostream& os, const distributionMap& map)
{
    os  << static_cast<const distributionMapBase&>(map) << token::NL
        << map.transformElements_ << token::NL
        << map.transformStart_;

    return os;
}


// ************************************************************************* //
