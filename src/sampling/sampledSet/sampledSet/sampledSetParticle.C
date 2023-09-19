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

#include "sampledSetParticle.H"
#include "sampledSetCloud.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::sampledSetParticle::sampledSetParticle
(
    const polyMesh& mesh,
    const point& position,
    const label celli,
    label& nLocateBoundaryHits,
    const label seti,
    const scalar setF,
    const scalar distance
)
:
    particle(mesh, position, celli, nLocateBoundaryHits),
    seti_(seti),
    setF_(setF),
    distance_(distance),
    havePosition0_(false),
    position0_(point::uniform(NaN))
{}


Foam::sampledSetParticle::sampledSetParticle(Istream& is, bool readFields)
:
    particle(is, readFields)
{
    if (readFields)
    {
        is  >> seti_ >> setF_ >> distance_ >> havePosition0_ >> position0_;
    }

    is.check(FUNCTION_NAME);
}


Foam::sampledSetParticle::sampledSetParticle
(
    const sampledSetParticle& p
)
:
    particle(p),
    seti_(p.seti_),
    setF_(p.setF_),
    distance_(p.distance_),
    havePosition0_(p.havePosition0_),
    position0_(p.position0_)
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::sampledSetParticle::store
(
    sampledSetCloud& cloud,
    trackingData& td
)
{
    td.positions_.append(position(td.mesh));
    td.distances_.append(distance_);
    td.cells_.append(cell());
    td.faces_.append(face());
}


void Foam::sampledSetParticle::storeFace
(
    sampledSetCloud& cloud,
    trackingData& td
)
{
    store(cloud, td);
}


void Foam::sampledSetParticle::storeCell
(
    sampledSetCloud& cloud,
    trackingData& td
)
{
    if (havePosition0_)
    {
        td.positions_.append((position(td.mesh) + position0_)/2);
        td.distances_.append(distance_ - mag(position(td.mesh) - position0_)/2);
        td.cells_.append(cell());
        td.faces_.append(-1);
    }

    havePosition0_ = true;
    position0_ = position(td.mesh);
}


bool Foam::sampledSetParticle::move
(
    sampledSetCloud& cloud,
    trackingData& td
)
{
    td.keepParticle = true;
    td.sendToProc = -1;

    while (td.keepParticle && td.sendToProc == -1 && seti_ < td.set_.size() - 1)
    {
        const vector s = td.set_[seti_ + 1] - td.set_[seti_];
        const scalar magS = mag(s);

        const scalar f = trackToFace(td.mesh, setF_*s, 0);
        distance_ += (1 - f)*setF_*magS;
        setF_ *= f;

        while (onInternalFace(td.mesh))
        {
            if (td.storeCells_) storeCell(cloud, td);
            if (td.storeFaces_ > 0) storeFace(cloud, td);

            hitFace(setF_*s, 0, cloud, td);

            const scalar f = trackToFace(td.mesh, setF_*s, 0);
            distance_ += (1 - f)*setF_*magS;
            setF_ *= f;
        }

        if (onFace())
        {
            if (td.storeCells_) storeCell(cloud, td);
            if (td.storeFaces_ > 0) storeFace(cloud, td);

            hitFace(setF_*s, 0, cloud, td);
        }

        if (onBoundaryFace(td.mesh) && setF_ < rootSmall)
        {
            if (td.storeSet_) store(cloud, td);
        }

        if (!onFace())
        {
            if (td.storeSet_) store(cloud, td);

            ++ seti_;
            setF_ = 1;
        }
    }

    return true;
}


void Foam::sampledSetParticle::hitWedgePatch
(
    sampledSetCloud&,
    trackingData& td
)
{
    seti_ = labelMax;
}


void Foam::sampledSetParticle::hitSymmetryPlanePatch
(
    sampledSetCloud&,
    trackingData& td
)
{
    seti_ = labelMax;
}


void Foam::sampledSetParticle::hitSymmetryPatch
(
    sampledSetCloud&,
    trackingData& td
)
{
    seti_ = labelMax;
}


void Foam::sampledSetParticle::hitCyclicPatch
(
    sampledSetCloud& cloud,
    trackingData& td
)
{
    seti_ = labelMax;
}


void Foam::sampledSetParticle::hitProcessorPatch
(
    sampledSetCloud& cloud,
    trackingData& td
)
{
    const processorPolyPatch& ppp =
        static_cast<const processorPolyPatch&>
        (
            td.mesh.boundaryMesh()[patch(td.mesh)]
        );

    if (ppp.transform().transformsPosition())
    {
        seti_ = labelMax;
    }
    else
    {
        particle::hitProcessorPatch(cloud, td);
    }
}


void Foam::sampledSetParticle::hitWallPatch
(
    sampledSetCloud&,
    trackingData& td
)
{
    seti_ = labelMax;
}


void Foam::sampledSetParticle::correctAfterParallelTransfer
(
    sampledSetCloud& cloud,
    trackingData& td
)
{
    particle::correctAfterParallelTransfer(cloud, td);

    if (td.storeFaces_ > 1) storeFace(cloud, td);
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<(Ostream& os, const sampledSetParticle& p)
{
    os  << static_cast<const particle&>(p)
        << token::SPACE << p.seti_
        << token::SPACE << p.setF_
        << token::SPACE << p.distance_
        << token::SPACE << p.havePosition0_
        << token::SPACE << p.position0_;

    os.check(FUNCTION_NAME);

    return os;
}


// ************************************************************************* //
