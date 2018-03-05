/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
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

\*---------------------------------------------------------------------------*/

#include "topoSet.H"
#include "mapPolyMesh.H"
#include "polyMesh.H"
#include "boundBox.H"
#include "Time.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(topoSet, 0);
    defineRunTimeSelectionTable(topoSet, word);
    defineRunTimeSelectionTable(topoSet, size);
    defineRunTimeSelectionTable(topoSet, set);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::topoSet> Foam::topoSet::New
(
    const word& setType,
    const polyMesh& mesh,
    const word& name,
    readOption r,
    writeOption w
)
{
    wordConstructorTable::iterator cstrIter =
        wordConstructorTablePtr_->find(setType);

    if (cstrIter == wordConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown set type " << setType
            << endl << endl
            << "Valid set types : " << endl
            << wordConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<topoSet>(cstrIter()(mesh, name, r, w));
}


Foam::autoPtr<Foam::topoSet> Foam::topoSet::New
(
    const word& setType,
    const polyMesh& mesh,
    const word& name,
    const label size,
    writeOption w
)
{
    sizeConstructorTable::iterator cstrIter =
        sizeConstructorTablePtr_->find(setType);

    if (cstrIter == sizeConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown set type " << setType
            << endl << endl
            << "Valid set types : " << endl
            << sizeConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<topoSet>(cstrIter()(mesh, name, size, w));
}


Foam::autoPtr<Foam::topoSet> Foam::topoSet::New
(
    const word& setType,
    const polyMesh& mesh,
    const word& name,
    const topoSet& set,
    writeOption w
)
{
    setConstructorTable::iterator cstrIter =
        setConstructorTablePtr_->find(setType);

    if (cstrIter == setConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown set type " << setType
            << endl << endl
            << "Valid set types : " << endl
            << setConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<topoSet>(cstrIter()(mesh, name, set, w));
}


Foam::fileName Foam::topoSet::localPath
(
    const polyMesh& mesh,
    const word& name
)
{
    return mesh.facesInstance()/mesh.dbDir()/polyMesh::meshSubDir/"sets"/name;
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

// Update stored cell numbers using map.
// Do in two passes to prevent allocation if nothing changed.
void Foam::topoSet::updateLabels(const labelList& map)
{
    // Iterate over map to see if anything changed
    bool changed = false;

    forAllConstIter(labelHashSet, *this, iter)
    {
        if ((iter.key() < 0) || (iter.key() > map.size()))
        {
            FatalErrorInFunction
                << "Illegal content " << iter.key() << " of set:" << name()
                << " of type " << type() << endl
                << "Value should be between 0 and " << map.size()-1
                << abort(FatalError);
        }

        const label newCelli = map[iter.key()];

        if (newCelli != iter.key())
        {
            changed = true;

            break;
        }
    }

    // Relabel (use second Map to prevent overlapping)
    if (changed)
    {
        labelHashSet newSet(2*size());

        forAllConstIter(labelHashSet, *this, iter)
        {
            const label newCelli = map[iter.key()];

            if (newCelli >= 0)
            {
                newSet.insert(newCelli);
            }
        }

        transfer(newSet);
    }
}


void Foam::topoSet::check(const label maxLabel)
{
    forAllConstIter(topoSet, *this, iter)
    {
        if ((iter.key() < 0) || (iter.key() > maxLabel))
        {
            FatalErrorInFunction
                << "Illegal content " << iter.key() << " of set:" << name()
                << " of type " << type() << endl
                << "Value should be between 0 and " << maxLabel
                << abort(FatalError);
        }
    }
}


// Write maxElem elements, starting at iter. Updates iter and elemI.
void Foam::topoSet::writeDebug
(
    Ostream& os,
    const label maxElem,
    topoSet::const_iterator& iter,
    label& elemI
) const
{
    label n = 0;

    for (; (iter != end()) && (n < maxElem); ++iter)
    {
        if ((n != 0) && ((n % 10) == 0))
        {
            os << endl;
        }
        os << iter.key() << ' ';

        n++;
        elemI++;
    }
}


// Write maxElem elements, starting at iter. Updates iter and elemI.
void Foam::topoSet::writeDebug
(
    Ostream& os,
    const pointField& coords,
    const label maxElem,
    topoSet::const_iterator& iter,
    label& elemI
) const
{
    label n = 0;

    for (; (iter != end()) && (n < maxElem); ++iter)
    {
        if ((n != 0) && ((n % 3) == 0))
        {
            os << endl;
        }
        os << iter.key() << coords[iter.key()] << ' ';

        n++;
        elemI++;
    }
}


void Foam::topoSet::writeDebug
(
    Ostream& os,
    const pointField& coords,
    const label maxLen
) const
{
    // Bounding box of contents.
    boundBox bb(pointField(coords, toc()), true);

    os  << "Set bounding box: min = "
        << bb.min() << "    max = " << bb.max() << " meters. " << endl << endl;

    label n = 0;

    topoSet::const_iterator iter = begin();

    if (size() <= maxLen)
    {
        writeDebug(os, coords, maxLen, iter, n);
    }
    else
    {
        label halfLen = maxLen/2;

        os  << "Size larger than " << maxLen << ". Printing first and last "
            << halfLen << " elements:" << endl << endl;

        writeDebug(os, coords, halfLen, iter, n);

        os  << nl << "  .." << nl << endl;

        for (; n < size() - halfLen; ++n)
        {
            ++iter;
        }

        writeDebug(os, coords, halfLen, iter, n);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::topoSet::topoSet(const IOobject& obj, const word& wantedType)
:
    regIOobject(obj)
{
    if
    (
        readOpt() == IOobject::MUST_READ
     || readOpt() == IOobject::MUST_READ_IF_MODIFIED
     || (
            readOpt() == IOobject::READ_IF_PRESENT
         && headerOk()
        )
    )
    {
        if (readStream(wantedType).good())
        {
            readStream(wantedType) >> static_cast<labelHashSet&>(*this);

            close();
        }
    }
}


Foam::topoSet::topoSet
(
    const polyMesh& mesh,
    const word& wantedType,
    const word& name,
    readOption r,
    writeOption w
)
:
    regIOobject
    (
        IOobject
        (
            name,
            mesh.time().findInstance
            (
                mesh.dbDir()/polyMesh::meshSubDir/"sets",
                word::null,
                r,  //IOobject::MUST_READ,
                mesh.facesInstance()
            ),
            polyMesh::meshSubDir/"sets",
            mesh,
            r,
            w
        )
    )
{
    if
    (
        readOpt() == IOobject::MUST_READ
     || readOpt() == IOobject::MUST_READ_IF_MODIFIED
     || (
            readOpt() == IOobject::READ_IF_PRESENT
         && headerOk()
        )
    )
    {
        if (readStream(wantedType).good())
        {
            readStream(wantedType) >> static_cast<labelHashSet&>(*this);

            close();
        }
    }
}


Foam::topoSet::topoSet
(
    const polyMesh& mesh,
    const word& name,
    const label size,
    writeOption w
)
:
    regIOobject
    (
        IOobject
        (
            name,
            mesh.time().findInstance
            (
                mesh.dbDir()/polyMesh::meshSubDir/"sets",
                word::null,
                IOobject::NO_READ,
                mesh.facesInstance()
            ),
            polyMesh::meshSubDir/"sets",
            mesh,
            IOobject::NO_READ,
            w
        )
    ),
    labelHashSet(size)
{}


Foam::topoSet::topoSet
(
    const polyMesh& mesh,
    const word& name,
    const labelHashSet& set,
    writeOption w
)
:
    regIOobject
    (
        IOobject
        (
            name,
            mesh.time().findInstance
            (
                mesh.dbDir()/polyMesh::meshSubDir/"sets",
                word::null,
                IOobject::NO_READ,
                mesh.facesInstance()
            ),
            polyMesh::meshSubDir/"sets",
            mesh,
            IOobject::NO_READ,
            w
        )
    ),
    labelHashSet(set)
{}


Foam::topoSet::topoSet(const IOobject& obj, const label size)
:
    regIOobject(obj),
    labelHashSet(size)
{}


Foam::topoSet::topoSet(const IOobject& obj, const labelHashSet& set)
:
    regIOobject(obj),
    labelHashSet(set)
{}



// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::topoSet::~topoSet()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::topoSet::invert(const label maxLen)
{
    // Keep copy of current set.
    labelHashSet currentSet(*this);

    clear();
    resize(2*(maxLen - currentSet.size()));

    for (label celli = 0; celli < maxLen; celli++)
    {
        if (!currentSet.found(celli))
        {
            insert(celli);
        }
    }
}


void Foam::topoSet::subset(const topoSet& set)
{
    // Keep copy of current set.
    labelHashSet currentSet(*this);

    clear();
    resize(2*min(currentSet.size(), set.size()));

    forAllConstIter(labelHashSet, currentSet, iter)
    {
        if (set.found(iter.key()))
        {
            // element present in both currentSet and set.
            insert(iter.key());
        }
    }
}


void Foam::topoSet::addSet(const topoSet& set)
{
    forAllConstIter(topoSet, set, iter)
    {
        insert(iter.key());
    }
}


void Foam::topoSet::deleteSet(const topoSet& set)
{
    forAllConstIter(topoSet, set, iter)
    {
        erase(iter.key());
    }
}


void Foam::topoSet::sync(const polyMesh&)
{
    NotImplemented;
}


void Foam::topoSet::writeDebug(Ostream& os, const label maxLen) const
{
    label n = 0;

    topoSet::const_iterator iter = begin();

    if (size() <= maxLen)
    {
        writeDebug(os, maxLen, iter, n);
    }
    else
    {
        label halfLen = maxLen/2;

        os  << "Size larger than " << maxLen << ". Printing first and last "
            << halfLen << " elements:" << endl << endl;

        writeDebug(os, halfLen, iter, n);

        os  << nl << "  .." << nl << endl;

        for (; n < size() - halfLen; ++n)
        {
            ++iter;
        }

        writeDebug(os, halfLen, iter, n);
    }
}


bool Foam::topoSet::writeData(Ostream& os) const
{
    return (os << *this).good();
}


void Foam::topoSet::updateMesh(const mapPolyMesh&)
{
    NotImplemented;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

void Foam::topoSet::operator=(const topoSet& rhs)
{
    labelHashSet::operator=(rhs);
}


// ************************************************************************* //
