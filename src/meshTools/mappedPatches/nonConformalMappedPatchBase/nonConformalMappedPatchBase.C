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

#include "nonConformalMappedPatchBase.H"
#include "nonConformalBoundary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(nonConformalMappedPatchBase, 0);
}


// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

bool Foam::nonConformalMappedPatchBase::calcOwner() const
{
    const polyPatch& pp = patch_.patch();

    const word& regionName = pp.boundaryMesh().mesh().name();
    const word& patchName = pp.name();

    if (regionName != nbrRegionName())
    {
        return regionName < nbrRegionName();
    }

    if (patchName != nbrPatchName())
    {
        return patchName < nbrPatchName();
    }

    FatalErrorInFunction
        << "Patch " << patch_.patch().name() << " of type "
        << nonConformalMappedPatchBase::typeName << " maps to itself"
        << exit(FatalError);

    return false;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::nonConformalMappedPatchBase::nonConformalMappedPatchBase
(
    const nonConformalPolyPatch& ncPp
)
:
    mappedPatchBaseBase(ncPp.patch()),
    patch_(ncPp),
    owner_(calcOwner()),
    intersectionIsValid_(0),
    intersection_(false)
{}


Foam::nonConformalMappedPatchBase::nonConformalMappedPatchBase
(
    const nonConformalPolyPatch& ncPp,
    const word& nbrRegionName,
    const word& nbrPatchName,
    const cyclicTransform& transform
)
:
    mappedPatchBaseBase(ncPp.patch(), nbrRegionName, nbrPatchName, transform),
    patch_(ncPp),
    owner_(calcOwner()),
    intersectionIsValid_(0),
    intersection_(false)
{}


Foam::nonConformalMappedPatchBase::nonConformalMappedPatchBase
(
    const nonConformalPolyPatch& ncPp,
    const dictionary& dict,
    const transformType tt
)
:
    mappedPatchBaseBase(ncPp.patch(), dict, tt),
    patch_(ncPp),
    owner_(calcOwner()),
    intersectionIsValid_(0),
    intersection_(false)
{}


Foam::nonConformalMappedPatchBase::nonConformalMappedPatchBase
(
    const nonConformalPolyPatch& ncPp,
    const nonConformalMappedPatchBase& mpb
)
:
    mappedPatchBaseBase(ncPp.patch(), mpb),
    patch_(ncPp),
    owner_(calcOwner()),
    intersectionIsValid_(0),
    intersection_(false)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::nonConformalMappedPatchBase::~nonConformalMappedPatchBase()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::nonConformalMappedPatchBase::owner() const
{
    return owner_;
}


const Foam::patchToPatches::intersection&
Foam::nonConformalMappedPatchBase::intersection() const
{
    if (!owner())
    {
        FatalErrorInFunction
            << "The non-conformal mapped intersection is only available to "
            << "the owner patch" << abort(FatalError);
    }

    const label bothIntersectionIsValid =
        min(intersectionIsValid_, nbrMappedPatch().intersectionIsValid_);

    const bool intersectionIsValid =
        (bothIntersectionIsValid == 2)
     || (
            bothIntersectionIsValid == 1
         && !mappedPatchBaseBase::moving
            (
                patch_.origPatch(),
                nbrMappedPatch().patch_.origPatch()
            )
        );

    if (!intersectionIsValid)
    {
        const polyMesh& mesh =
            mappedPatchBaseBase::patch_.boundaryMesh().mesh();

        const nonConformalBoundary& ncb = nonConformalBoundary::New(mesh);

        intersection_.update
        (
            patch_.origPatch(),
            ncb.patchPointNormals(patch_.origPatchIndex()),
            nbrMappedPatch().patch_.origPatch(),
            transform_.transform()
        );

        intersectionIsValid_ = 2;

        nbrMappedPatch().intersectionIsValid_ = 2;
    }

    return intersection_;
}


void Foam::nonConformalMappedPatchBase::clearOut(const bool move)
{
    if (move && moveUpdate_ == moveUpdate::never)
    {
        // Do nothing
    }
    else if (move && moveUpdate_ == moveUpdate::detect)
    {
        intersectionIsValid_ = min(intersectionIsValid_, 1);
    }
    else
    {
        intersectionIsValid_ = 0;
    }
}


void Foam::nonConformalMappedPatchBase::write(Ostream& os) const
{
    mappedPatchBaseBase::write(os);
}


// ************************************************************************* //
