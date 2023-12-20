/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2021-2023 OpenFOAM Foundation
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

#include "nonConformalPolyPatch.H"
#include "polyMesh.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(nonConformalPolyPatch, 0);
}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

void Foam::nonConformalPolyPatch::validateSize() const
{
    if (patch_.size() != 0)
    {
        FatalErrorInFunction
            << "Patch " << patch_.name() << " has " << patch_.size()
            << " faces. Patches of type " << patch_.type()
            << " must have zero faces." << exit(FatalError);
    }
}


// * * * * * * * * * * * *  Protected Member Functions * * * * * * * * * * * //

void Foam::nonConformalPolyPatch::rename(const wordList& newNames)
{
    if (origPatchIndex_ != -1)
    {
        origPatchName_ = newNames[origPatchIndex_];
    }
    else
    {
        FatalErrorInFunction
            << "Cannot rename " << nonConformalPolyPatch::typeName
            << " without the original patch index"
            << exit(FatalError);
    }
}


void Foam::nonConformalPolyPatch::reorder(const labelUList& newToOldIndex)
{
    if (origPatchIndex_ != -1)
    {
        origPatchIndex_ = findIndex(newToOldIndex, origPatchIndex_);
    }
}


// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //

Foam::nonConformalPolyPatch::nonConformalPolyPatch(const polyPatch& patch)
:
    patch_(patch),
    origPatchName_(word::null),
    origPatchIndex_(-1)
{
    validateSize();
}


Foam::nonConformalPolyPatch::nonConformalPolyPatch
(
    const polyPatch& patch,
    const word& origPatchName
)
:
    patch_(patch),
    origPatchName_(origPatchName),
    origPatchIndex_(-1)
{
    validateSize();
}


Foam::nonConformalPolyPatch::nonConformalPolyPatch
(
    const polyPatch& patch,
    const dictionary& dict
)
:
    patch_(patch),
    origPatchName_(dict.lookup<word>("originalPatch")),
    origPatchIndex_(-1)
{
    validateSize();
}


Foam::nonConformalPolyPatch::nonConformalPolyPatch
(
    const polyPatch& patch,
    const nonConformalPolyPatch& ncPatch
)
:
    patch_(patch),
    origPatchName_(ncPatch.origPatchName_),
    origPatchIndex_(-1)
{
    validateSize();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::nonConformalPolyPatch::~nonConformalPolyPatch()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::polyPatch& Foam::nonConformalPolyPatch::patch() const
{
    return patch_;
}


const Foam::word& Foam::nonConformalPolyPatch::origPatchName() const
{
    return origPatchName_;
}


Foam::label Foam::nonConformalPolyPatch::origPatchIndex() const
{
    if (origPatchIndex_ == -1)
    {
        origPatchIndex_ = patch_.boundaryMesh().findIndex(origPatchName());

        if (origPatchIndex_ == -1)
        {
            FatalErrorInFunction
                << "Illegal neighbourPatch name " << origPatchName()
                << endl << "Valid patch names are "
                << patch_.boundaryMesh().names()
                << exit(FatalError);
        }

        const polyPatch& p = patch_.boundaryMesh()[origPatchIndex_];

        if (isA<nonConformalPolyPatch>(p))
        {
            FatalErrorInFunction
                << "The originalPatch for the "
                << patch_.type() << " patch " << patch_.name() << " is "
                << p.name() << " which is also of "
                << nonConformalPolyPatch::typeName << " type. This is not "
                << "allowed. The originalPatch must be of a non-"
                << nonConformalPolyPatch::typeName << " type."
                << exit(FatalError);
        }
    }

    return origPatchIndex_;
}


const Foam::polyPatch& Foam::nonConformalPolyPatch::origPatch() const
{
    return patch_.boundaryMesh()[origPatchIndex()];
}


void Foam::nonConformalPolyPatch::write(Ostream& os) const
{
    writeEntry(os, "originalPatch", origPatchName_);
}


// ************************************************************************* //
