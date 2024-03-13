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

#include "zone.H"
#include "IOstream.H"
#include "demandDrivenData.H"
#include "HashSet.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(zone, 0);
}


// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

const Foam::Map<Foam::label>& Foam::zone::lookupMap() const
{
    if (!lookupMapPtr_)
    {
        calcLookupMap();
    }

    return *lookupMapPtr_;
}


void Foam::zone::calcLookupMap() const
{
    if (debug)
    {
        InfoInFunction << "Calculating lookup map" << endl;
    }

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

    if (debug)
    {
        InfoInFunction << "Finished calculating lookup map" << endl;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::zone::zone
(
    const word& name,
    const labelUList& indices
)
:
    labelList(indices),
    name_(name),
    lookupMapPtr_(nullptr)
{}


Foam::zone::zone
(
    const word& name,
    labelList&& indices
)
:
    labelList(move(indices)),
    name_(name),
    lookupMapPtr_(nullptr)
{}


Foam::zone::zone
(
    const word& name,
    const dictionary& dict,
    const word& labelsName
)
:
    labelList(dict.lookup(labelsName)),
    name_(name),
    lookupMapPtr_(nullptr)
{}


Foam::zone::zone
(
    const zone& z,
    const labelUList& indices
)
:
    labelList(indices),
    name_(z.name()),
    lookupMapPtr_(nullptr)
{}


Foam::zone::zone
(
    const zone& z,
    labelList&& indices
)
:
    labelList(move(indices)),
    name_(z.name()),
    lookupMapPtr_(nullptr)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::zone::~zone()
{
    clearAddressing();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::label Foam::zone::localIndex(const label globalIndex) const
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


void Foam::zone::clearAddressing()
{
    deleteDemandDrivenData(lookupMapPtr_);
}


bool Foam::zone::checkDefinition(const label maxSize, const bool report) const
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


void Foam::zone::insert(const labelHashSet& newIndices)
{
    labelHashSet indices(*this);
    indices.insert(newIndices);
    labelList::operator=(indices.sortedToc());
}


void Foam::zone::mapMesh(const polyMeshMap&)
{
    clearAddressing();
}


void Foam::zone::distribute(const polyDistributionMap&)
{
    clearAddressing();
}


void Foam::zone::write(Ostream& os) const
{
    os  << nl << name_
        << nl << static_cast<const labelList&>(*this);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

void Foam::zone::operator=(const zone& zn)
{
    clearAddressing();
    labelList::operator=(zn);
}


void Foam::zone::operator=(zone&& zn)
{
    clearAddressing();
    labelList::operator=(move(zn));
}


void Foam::zone::operator=(const labelUList& indices)
{
    clearAddressing();
    labelList::operator=(indices);
}


void Foam::zone::operator=(labelList&& indices)
{
    clearAddressing();
    labelList::operator=(move(indices));
}


// * * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<(Ostream& os, const zone& z)
{
    z.write(os);
    os.check("Ostream& operator<<(Ostream& f, const zone& z");
    return os;
}


// ************************************************************************* //
