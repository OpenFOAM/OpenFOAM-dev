/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2024 OpenFOAM Foundation
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

#include "Zone.H"
#include "HashSet.H"
#include "pointField.H"
#include "demandDrivenData.H"

// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

template<class ZoneType, class ZonesType>
template<class Type>
void Foam::Zone<ZoneType, ZonesType>::select(const Type& zone)
{
    const pointField& ctrs = zone.meshCentres();

    labelList& indices = *this;
    indices.setSize(ctrs.size());

    label nInZone = 0;
    forAll(ctrs, i)
    {
        if (zone.contains(ctrs[i]))
        {
            indices[nInZone++] = i;
        }
    }

    indices.setSize(nInZone);
}


template<class ZoneType, class ZonesType>
void Foam::Zone<ZoneType, ZonesType>::calcLookupMap() const
{
    if (lookupMapPtr_)
    {
        FatalErrorInFunction
            << "Lookup map already calculated" << nl
            << abort(FatalError);
    }

    const labelList& indices = *this;

    lookupMapPtr_ = new Map<label>(2*indices.size());
    Map<label>& lm = *lookupMapPtr_;

    forAll(indices, i)
    {
        lm.insert(indices[i], i);
    }
}


template<class ZoneType, class ZonesType>
void Foam::Zone<ZoneType, ZonesType>::topoChange
(
    const labelList& map,
    const labelList& reverseMap
)
{
    clearAddressing();

    labelHashSet indices;

    forAll(map, i)
    {
        if (map[i] >= 0 && localIndex(map[i]) != -1)
        {
            indices.insert(i);
        }
    }

    forAll(reverseMap, i)
    {
        if (reverseMap[i] >= 0 && localIndex(i) != -1)
        {
            indices.insert(reverseMap[i]);
        }
    }

    labelList::operator=(indices.sortedToc());
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ZoneType, class ZonesType>
Foam::Zone<ZoneType, ZonesType>::Zone
(
    const word& name,
    const labelUList& indices,
    const ZonesType& zones
)
:
    labelList(indices),
    name_(name),
    zones_(zones),
    lookupMapPtr_(nullptr)
{}


template<class ZoneType, class ZonesType>
Foam::Zone<ZoneType, ZonesType>::Zone
(
    const word& name,
    labelList&& indices,
    const ZonesType& zones
)
:
    labelList(move(indices)),
    name_(name),
    zones_(zones),
    lookupMapPtr_(nullptr)
{}


template<class ZoneType, class ZonesType>
Foam::Zone<ZoneType, ZonesType>::Zone
(
    const word& name,
    const dictionary& dict,
    const ZonesType& zones
)
:
    labelList(dict.lookup(ZoneType::labelsName)),
    name_(name),
    zones_(zones),
    lookupMapPtr_(nullptr)
{}


template<class ZoneType, class ZonesType>
Foam::Zone<ZoneType, ZonesType>::Zone
(
    const Zone& z,
    const labelUList& indices,
    const ZonesType& zones
)
:
    labelList(indices),
    name_(z.name()),
    zones_(zones),
    lookupMapPtr_(nullptr)
{}


template<class ZoneType, class ZonesType>
Foam::Zone<ZoneType, ZonesType>::Zone
(
    const Zone& z,
    labelList&& indices,
    const ZonesType& zones
)
:
    labelList(move(indices)),
    name_(z.name()),
    zones_(zones),
    lookupMapPtr_(nullptr)
{}


template<class ZoneType, class ZonesType>
Foam::autoPtr<ZoneType> Foam::Zone<ZoneType, ZonesType>::New
(
    const word& name,
    const dictionary& dict,
    const ZonesType& mz
)
{
    if (ZoneType::debug)
    {
        InfoInFunction
            << "Constructing " << ZoneType::typeName << " " << name << endl;
    }

    const word type(dict.lookup("type"));

    typename dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(type);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalIOErrorInFunction
        (
            dict
        )   << "Unknown " << ZoneType::typeName << " type "
            << type << nl << nl
            << "Valid " << ZoneType::typeName << " types are:" << nl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalIOError);
    }

    return autoPtr<ZoneType>(cstrIter()(name, dict, mz));
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class ZoneType, class ZonesType>
Foam::Zone<ZoneType, ZonesType>::~Zone()
{
    clearAddressing();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ZoneType, class ZonesType>
const ZonesType& Foam::Zone<ZoneType, ZonesType>::zones() const
{
    return zones_;
}


template<class ZoneType, class ZonesType>
Foam::label Foam::Zone<ZoneType, ZonesType>::localIndex
(
    const label globalIndex
) const
{
    const Map<label>& lm = lookupMap();

    Map<label>::const_iterator lmIter = lm.find(globalIndex);

    if (lmIter == lm.end())
    {
        return -1;
    }
    else
    {
        return lmIter();
    }
}


template<class ZoneType, class ZonesType>
const Foam::Map<Foam::label>& Foam::Zone<ZoneType, ZonesType>::lookupMap() const
{
    if (!lookupMapPtr_)
    {
        calcLookupMap();
    }

    return *lookupMapPtr_;
}


template<class ZoneType, class ZonesType>
void Foam::Zone<ZoneType, ZonesType>::clearAddressing()
{
    deleteDemandDrivenData(lookupMapPtr_);
}


template<class ZoneType, class ZonesType>
bool Foam::Zone<ZoneType, ZonesType>::checkDefinition
(
    const label maxSize,
    const bool report
) const
{
    const labelList& indices = *this;

    bool hasError = false;

    // To check for duplicate entries
    labelHashSet elems(size());

    forAll(indices, i)
    {
        if (indices[i] < 0 || indices[i] >= maxSize)
        {
            hasError = true;

            if (report)
            {
                SeriousErrorInFunction
                    << "Zone " << name_
                    << " contains invalid index label " << indices[i] << nl
                    << "Valid index labels are 0.."
                    << maxSize-1 << endl;
            }
            else
            {
                // w/o report - can stop checking now
                break;
            }
        }
        else if (!elems.insert(indices[i]))
        {
            if (report)
            {
                WarningInFunction
                    << "Zone " << name_
                    << " contains duplicate index label " << indices[i] << endl;
            }
        }
    }

    return hasError;
}


template<class ZoneType, class ZonesType>
void Foam::Zone<ZoneType, ZonesType>::insert(const labelHashSet& newIndices)
{
    labelHashSet indices(*this);
    indices.insert(newIndices);
    labelList::operator=(indices.sortedToc());
}


template<class ZoneType, class ZonesType>
void Foam::Zone<ZoneType, ZonesType>::swap(Zone& z)
{
    clearAddressing();
    labelList::swap(z);
}


template<class ZoneType, class ZonesType>
void Foam::Zone<ZoneType, ZonesType>::mapMesh(const polyMeshMap&)
{
    clearAddressing();
}


template<class ZoneType, class ZonesType>
void Foam::Zone<ZoneType, ZonesType>::distribute(const polyDistributionMap&)
{
    clearAddressing();
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class ZoneType, class ZonesType>
void Foam::Zone<ZoneType, ZonesType>::operator=(const Zone& zn)
{
    clearAddressing();
    labelList::operator=(zn);
}


template<class ZoneType, class ZonesType>
void Foam::Zone<ZoneType, ZonesType>::operator=(Zone&& zn)
{
    clearAddressing();
    labelList::operator=(move(zn));
}


template<class ZoneType, class ZonesType>
void Foam::Zone<ZoneType, ZonesType>::operator=(const labelUList& indices)
{
    clearAddressing();
    labelList::operator=(indices);
}


template<class ZoneType, class ZonesType>
void Foam::Zone<ZoneType, ZonesType>::operator=(labelList&& indices)
{
    clearAddressing();
    labelList::operator=(move(indices));
}


// * * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * //

template<class ZoneType, class ZonesType>
Foam::Ostream& Foam::operator<<(Ostream& os, const Zone<ZoneType, ZonesType>& z)
{
    z.writeDict(os);
    os.check("Ostream& operator<<(Ostream& f, const Zone& z");
    return os;
}


// ************************************************************************* //
