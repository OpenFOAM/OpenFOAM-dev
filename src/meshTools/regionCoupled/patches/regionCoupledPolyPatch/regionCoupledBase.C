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

#include "regionCoupledBase.H"
#include "SubField.H"
#include "polyMesh.H"
#include "Time.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(regionCoupledBase, 0);
}


// * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * * //

void Foam::regionCoupledBase::resetAMI() const
{
    if (owner())
    {
        AMIPtr_.clear();

        const polyPatch& nbr = refCast<const polyPatch>(neighbPatch());
        pointField nbrPoints = nbr.localPoints();

        if (debug)
        {
            const Time& t = patch_.boundaryMesh().mesh().time();
            OFstream os(t.path()/patch_.name() + "_neighbourPatch-org.obj");
            meshTools::writeOBJ(os, nbr.localFaces(), nbrPoints);
        }

        // transform neighbour patch to local system
        // transformPosition(nbrPoints);
        primitivePatch nbrPatch0
        (
            SubList<face>
            (
                nbr.localFaces(),
                nbr.size()
            ),
            nbrPoints
        );

        if (debug)
        {
            const Time& t = patch_.boundaryMesh().mesh().time();
            OFstream osN(t.path()/patch_.name() + "_neighbourPatch-trans.obj");
            meshTools::writeOBJ(osN, nbrPatch0.localFaces(), nbrPoints);

            OFstream osO(t.path()/patch_.name() + "_ownerPatch.obj");
            meshTools::writeOBJ
            (
                osO,
                patch_.localFaces(),
                patch_.localPoints()
            );
        }

        // Construct/apply AMI interpolation to determine addressing and weights
        AMIPtr_.reset
        (
            new AMIPatchToPatchInterpolation
            (
                patch_,
                nbrPatch0,
                surfPtr(),
                faceAreaIntersect::tmMesh,
                true,
                AMIPatchToPatchInterpolation::imFaceAreaWeight,
                -1,
                AMIReverse_
            )
        );

        if (debug)
        {
            Pout<< "regionCoupledBase : " << patch_.name()
                << " constructed AMI with " << nl
                << "    " << ":srcAddress:" << AMIPtr_().srcAddress().size()
                << nl
                << "    " << " tgAddress :" << AMIPtr_().tgtAddress().size()
                << nl << endl;
        }
    }
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::regionCoupledBase::clearGeom()
{
    AMIPtr_.clear();
    surfPtr_.clear();
}


// * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * * //

Foam::regionCoupledBase::regionCoupledBase
(
    const polyPatch& pp
)
:
    patch_(pp),
    nbrPatchName_(word::null),
    nbrPatchID_(-1),
    nbrRegionName_(word::null),
    sameRegion_(false),
    AMIPtr_(nullptr),
    AMIReverse_(false),
    surfPtr_(nullptr),
    surfDict_(fileName("surface"))
{}


Foam::regionCoupledBase::regionCoupledBase
(
    const polyPatch& pp,
    const dictionary& dict
)
:
    patch_(pp),
    nbrPatchName_(dict.lookup("neighbourPatch")),
    nbrPatchID_(-1),
    nbrRegionName_(dict.lookup("neighbourRegion")),
    sameRegion_(nbrRegionName_ == patch_.boundaryMesh().mesh().name()),
    AMIPtr_(nullptr),
    AMIReverse_(dict.lookupOrDefault<bool>("flipNormals", false)),
    surfPtr_(nullptr),
    surfDict_(dict.subOrEmptyDict("surface"))
{}


Foam::regionCoupledBase::regionCoupledBase
(
    const polyPatch& pp,
    const regionCoupledBase& mpb
)
:
    patch_(pp),
    nbrPatchName_(mpb.nbrPatchName_),
    nbrPatchID_(mpb.nbrPatchID_),
    nbrRegionName_(mpb.nbrRegionName_),
    sameRegion_(mpb.sameRegion_),
    AMIPtr_(nullptr),
    AMIReverse_(mpb.AMIReverse_),
    surfPtr_(mpb.surfPtr_),
    surfDict_(mpb.surfDict_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::regionCoupledBase::~regionCoupledBase()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::label Foam::regionCoupledBase::neighbPatchID() const
{
    if (nbrPatchID_ == -1)
    {
        if
        (
            patch_.boundaryMesh().mesh().time().foundObject<polyMesh>
            (
                nbrRegionName_
            )
        )
        {
            const polyMesh& mesh =
                patch_.boundaryMesh().mesh().time().lookupObject<polyMesh>
                (
                    nbrRegionName_
                );

            nbrPatchID_ = mesh.boundaryMesh().findPatchID(nbrPatchName_);

            if (nbrPatchID_ == -1)
            {
                FatalErrorInFunction
                    << "Illegal neighbourPatch name " << nbrPatchName_
                    << nl << "Valid patch names are "
                    << mesh.boundaryMesh().names()
                    << exit(FatalError);
            }

            // Check that it is a cyclic AMI patch
            const regionCoupledBase& nbrPatch =
                refCast<const regionCoupledBase>
                (
                    mesh.boundaryMesh()[nbrPatchID_]
                );

            if (nbrPatch.nbrPatchName() != patch_.name())
            {
                WarningInFunction
                    << "Patch " << patch_.name()
                    << " specifies neighbour patch " << nbrPatchName()
                    << nl << " but that in return specifies "
                    << nbrPatch.nbrPatchName() << endl;
            }
        }
    }

    return nbrPatchID_;
}


bool Foam::regionCoupledBase::owner() const
{
    if (nbrRegionName_ == patch_.boundaryMesh().mesh().name())
    {
        return patch_.index() < neighbPatchID();
    }
    else
    {
        return patch_.boundaryMesh().mesh().name() < nbrRegionName_;
    }
}


const Foam::autoPtr<Foam::searchableSurface>&Foam::regionCoupledBase::
surfPtr() const
{
    const word surfType(surfDict_.lookupOrDefault<word>("type", "none"));

    if (!surfPtr_.valid() && owner() && surfType != "none")
    {
        word surfName(surfDict_.lookupOrDefault("name", patch_.name()));

        const polyMesh& mesh = patch_.boundaryMesh().mesh();

        surfPtr_ =
            searchableSurface::New
            (
                surfType,
                IOobject
                (
                    surfName,
                    mesh.time().constant(),
                    "triSurface",
                    mesh,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                ),
                surfDict_
            );
    }

    return surfPtr_;
}


const Foam::AMIPatchToPatchInterpolation& Foam::regionCoupledBase::AMI() const
{
    if (!owner())
    {
        FatalErrorInFunction
            << "AMI interpolator only available to owner patch"
            << abort(FatalError);
    }

    if (!AMIPtr_.valid())
    {
        resetAMI();
    }

    return AMIPtr_();
}


const Foam::regionCoupledBase&
Foam::regionCoupledBase::neighbPatch() const
{
    const polyMesh& mesh =
        patch_.boundaryMesh().mesh().time().lookupObject<polyMesh>
        (
            nbrRegionName_
        );

    const polyPatch& pp = mesh.boundaryMesh()[neighbPatchID()];
    return refCast<const regionCoupledBase>(pp);
}


bool Foam::regionCoupledBase::order
(
    PstreamBuffers& pBufs,
    const primitivePatch& pp,
    labelList& faceMap,
    labelList& rotation
) const
{
    faceMap.setSize(pp.size());
    faceMap = -1;

    rotation.setSize(pp.size());
    rotation = 0;

    return false;
}


void Foam::regionCoupledBase::write(Ostream& os) const
{
    os.writeKeyword("neighbourPatch") << nbrPatchName_
    << token::END_STATEMENT << nl;
    os.writeKeyword("neighbourRegion") << nbrRegionName_
    << token::END_STATEMENT << nl;

    if (AMIReverse_)
    {
        os.writeKeyword("flipNormals") << AMIReverse_
            << token::END_STATEMENT << nl;
    }

    if (!surfDict_.empty())
    {
        os.writeKeyword(surfDict_.dictName());
        os  << surfDict_;
    }
}


// ************************************************************************* //
