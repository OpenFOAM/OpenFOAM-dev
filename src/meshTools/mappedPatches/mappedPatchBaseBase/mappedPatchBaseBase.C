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

#include "mappedPatchBaseBase.H"
#include "Time.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(mappedPatchBaseBase, 0);
}


namespace Foam
{
    template<>
    const char*
        NamedEnum<mappedPatchBaseBase::moveUpdate, 3>::names[] =
        {"always", "detect", "never"};
}

const Foam::NamedEnum<Foam::mappedPatchBaseBase::moveUpdate, 3>
    Foam::mappedPatchBaseBase::moveUpdateNames_;


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::mappedPatchBaseBase::mappedPatchBaseBase(const polyPatch& pp)
:
    patch_(pp),
    coupleGroup_(),
    nbrRegionName_(patch_.boundaryMesh().mesh().name()),
    nbrPatchName_(patch_.name()),
    transform_(true),
    moveUpdate_(moveUpdate::always)
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
    transform_(transform),
    moveUpdate_(moveUpdate::always)
{}


Foam::mappedPatchBaseBase::mappedPatchBaseBase
(
    const polyPatch& pp,
    const dictionary& dict,
    const transformType tt
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
        tt == transformType::none ? cyclicTransform(true)
      : tt == transformType::defaultNone ? cyclicTransform(dict, true)
      : tt == transformType::specified ? cyclicTransform(dict, false)
      : cyclicTransform()
    ),
    moveUpdate_
    (
        dict.found("moveUpdate")
      ? moveUpdateNames_.read(dict.lookup("moveUpdate"))
      : dict.found("reMapAfterMove") // <-- backwards compatibility
      ? (
            dict.lookup<bool>("reMapAfterMove")
          ? moveUpdate::always
          : moveUpdate::never
        )
      : moveUpdate::always
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

    if (tt == transformType::none)
    {
        const cyclicTransform::transformTypes cttt =
            cyclicTransform(dict, true).transformType();

        if (cttt != cyclicTransform::NONE)
        {
            FatalIOErrorInFunction(dict)
                << word(cyclicTransform::transformTypeNames[cttt]).capitalise()
                << " transform specified for patch '" << patch_.name()
                << "' in region '" << patch_.boundaryMesh().mesh().name()
                << "'. This patch does not support transformed mapping."
                << exit(FatalIOError);
        }
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
    transform_(mpb.transform_),
    moveUpdate_(mpb.moveUpdate_)
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


bool Foam::mappedPatchBaseBase::moving
(
    const polyPatch& patch,
    const polyPatch& nbrPatch
)
{
    const polyMesh& mesh = patch.boundaryMesh().mesh();
    const polyMesh& nbrMesh = nbrPatch.boundaryMesh().mesh();

    if (!mesh.moving() && !nbrMesh.moving()) return false;

    auto localPatchMoving = [](const polyPatch& patch)
    {
        const polyMesh& mesh = patch.boundaryMesh().mesh();

        if (!mesh.moving()) return false;

        const pointField& points = mesh.points();
        const pointField& oldPoints = mesh.oldPoints();

        // !!! Note that this is an exact binary comparison of the point
        // positions. This isn't typical, but in most cases where this is
        // relevant the mesh motion is being switched off entirely for a period
        // of time, or the motion is only present in another part of the mesh,
        // and the points in question therefore do not change their value in
        // any way. This prevents the need for a more substantial calculation
        // procedure and the specification of tolerances and similar.

        return
            UIndirectList<point>(points, patch.meshPoints())
         != UIndirectList<point>(oldPoints, patch.meshPoints());
    };

    return
        returnReduce
        (
            localPatchMoving(patch) || localPatchMoving(nbrPatch),
            orOp<bool>()
        );
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

    writeEntryIfDifferent<word>
    (
        os,
        "moveUpdate",
        moveUpdateNames_[moveUpdate::always],
        moveUpdateNames_[moveUpdate_]
    );
}


// ************************************************************************* //
