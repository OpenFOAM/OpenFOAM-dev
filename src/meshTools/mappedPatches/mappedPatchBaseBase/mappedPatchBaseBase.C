/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2023 OpenFOAM Foundation
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

#include "mappedPatchBaseBase.H"
#include "Time.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(mappedPatchBaseBase, 0);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::mappedPatchBaseBase::mappedPatchBaseBase(const polyPatch& pp)
:
    patch_(pp),
    coupleGroup_(),
    nbrRegionName_(patch_.boundaryMesh().mesh().name()),
    nbrPatchName_(patch_.name()),
    transform_(true)
{}


Foam::mappedPatchBaseBase::mappedPatchBaseBase
(
    const polyPatch& pp,
    const word& nbrRegionName,
    const word& nbrPatchName,
    const cyclicTransform& transform
)
:
    patch_(pp),
    coupleGroup_(),
    nbrRegionName_(nbrRegionName),
    nbrPatchName_(nbrPatchName),
    transform_(transform)
{}


Foam::mappedPatchBaseBase::mappedPatchBaseBase
(
    const polyPatch& pp,
    const dictionary& dict,
    const bool transformIsNone
)
:
    patch_(pp),
    coupleGroup_(dict),
    nbrRegionName_
    (
        coupleGroup_.valid() ? word::null
      : dict.lookupOrDefaultBackwardsCompatible<word>
        (
            {"neighbourRegion", "sampleRegion"},
            pp.boundaryMesh().mesh().name()
        )
    ),
    nbrPatchName_
    (
        coupleGroup_.valid() ? word::null
      : dict.lookupOrDefault<bool>("samePatch", false) ? pp.name()
      : dict.lookupBackwardsCompatible<word>({"neighbourPatch", "samplePatch"})
    ),
    transform_
    (
        transformIsNone
      ? cyclicTransform(true)
      : cyclicTransform(dict, false)
    )
{
    const bool haveCoupleGroup = coupleGroup_.valid();

    const bool haveNbrRegion =
        dict.found("neighbourRegion") || dict.found("sampleRegion");
    const bool haveNbrPatch =
        dict.found("neighbourPatch") || dict.found("samplePatch");

    const bool isSamePatch = dict.lookupOrDefault<bool>("samePatch", false);

    if ((haveNbrRegion || haveNbrPatch || isSamePatch) && haveCoupleGroup)
    {
        FatalIOErrorInFunction(dict)
            << "Either neighbourRegion/Patch information or a coupleGroup "
            << "should be specified, not both" << exit(FatalIOError);
    }

    if (haveNbrPatch && isSamePatch)
    {
        FatalIOErrorInFunction(dict)
            << "Either a neighbourPatch should be specified, or samePatch "
            << "should be set to true, not both" << exit(FatalIOError);
    }
}


Foam::mappedPatchBaseBase::mappedPatchBaseBase
(
    const polyPatch& pp,
    const mappedPatchBaseBase& mpb
)
:
    patch_(pp),
    coupleGroup_(mpb.coupleGroup_),
    nbrRegionName_(mpb.nbrRegionName_),
    nbrPatchName_(mpb.nbrPatchName_),
    transform_(mpb.transform_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::mappedPatchBaseBase::~mappedPatchBaseBase()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::mappedPatchBaseBase::haveNbr() const
{
    const polyMesh& mesh = patch_.boundaryMesh().mesh();

    return mesh.time().foundObject<polyMesh>(nbrRegionName());
}


const Foam::polyMesh& Foam::mappedPatchBaseBase::nbrMesh() const
{
    const polyMesh& mesh = patch_.boundaryMesh().mesh();

    return mesh.time().lookupObject<polyMesh>(nbrRegionName());
}


const Foam::polyPatch& Foam::mappedPatchBaseBase::nbrPolyPatch() const
{
    const polyMesh& nbrMesh = this->nbrMesh();

    const label patchi = nbrMesh.boundaryMesh().findIndex(nbrPatchName());

    if (patchi == -1)
    {
        FatalErrorInFunction
            << "Cannot find patch " << nbrPatchName()
            << " in region " << nbrRegionName() << endl
            << "Valid patches are " << nbrMesh.boundaryMesh().names()
            << exit(FatalError);
    }

    return nbrMesh.boundaryMesh()[patchi];
}


bool Foam::mappedPatchBaseBase::specified(const dictionary& dict)
{
    return
        dict.found("coupleGroup")
     || dict.found("neighbourRegion")
     || dict.found("sampleRegion")
     || dict.found("neighbourPatch")
     || dict.found("samplePatch")
     || dict.found("samePatch");
}


void Foam::mappedPatchBaseBase::write(Ostream& os) const
{
    writeEntryIfDifferent(os, "neighbourRegion", word::null, nbrRegionName_);
    writeEntryIfDifferent(os, "neighbourPatch", word::null, nbrPatchName_);

    coupleGroup_.write(os);

    transform_.write(os);
}


// ************************************************************************* //
