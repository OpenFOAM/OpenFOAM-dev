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

#include "trackedParticle.H"


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::trackedParticle::trackedParticle
(
    const polyMesh& mesh,
    const vector& position,
    const label celli,
    const point& end,
    const label level,
    const label i,
    const label j,
    const label k
)
:
    particle(mesh, position, celli),
    start_(this->position(mesh)),
    end_(end),
    level_(level),
    i_(i),
    j_(j),
    k_(k)
{}


Foam::trackedParticle::trackedParticle(Istream& is, bool readFields)
:
    particle(is, readFields)
{
    if (readFields)
    {
        if (is.format() == IOstream::ASCII)
        {
            is >> start_ >> end_;
            level_ = readLabel(is);
            i_ = readLabel(is);
            j_ = readLabel(is);
            k_ = readLabel(is);
        }
        else
        {
            is.read
            (
                reinterpret_cast<char*>(&start_),
                sizeof(start_) + sizeof(end_) + sizeof(level_)
              + sizeof(i_) + sizeof(j_) + sizeof(k_)
            );
        }
    }

    // Check state of Istream
    is.check
    (
        "trackedParticle::trackedParticle"
        "(const Cloud<trackedParticle>&, Istream&, bool)"
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::trackedParticle::move
(
    Cloud<trackedParticle>& cloud,
    trackingData& td
)
{
    const scalar tEnd = (1.0 - stepFraction())*td.maxTrackLen_;

    if (tEnd <= small && onBoundaryFace(td.mesh))
    {
        // This is a hack to handle particles reaching their endpoint on a
        // processor boundary. This prevents a particle being endlessly
        // transferred backwards and forwards across the interface.
        td.keepParticle = false;
        td.sendToProc = -1;
    }
    else
    {
        td.keepParticle = true;
        td.sendToProc = -1;

        while (td.keepParticle && td.sendToProc == -1 && stepFraction() < 1)
        {
            // mark visited cell with max level.
            td.maxLevel_[cell()] = max(td.maxLevel_[cell()], level_);

            const scalar f = 1 - stepFraction();
            const vector s = end_ - start_;
            trackToAndHitFace(f*s, f, cloud, td);
        }
    }

    return td.keepParticle;
}


void Foam::trackedParticle::hitWedgePatch
(
    Cloud<trackedParticle>&,
    trackingData& td
)
{
    // Remove particle
    td.keepParticle = false;
}


void Foam::trackedParticle::hitSymmetryPlanePatch
(
    Cloud<trackedParticle>&,
    trackingData& td
)
{
    // Remove particle
    td.keepParticle = false;
}


void Foam::trackedParticle::hitSymmetryPatch
(
    Cloud<trackedParticle>&,
    trackingData& td
)
{
    // Remove particle
    td.keepParticle = false;
}


void Foam::trackedParticle::hitCyclicPatch
(
    Cloud<trackedParticle>&,
    trackingData& td
)
{
    // Remove particle
    td.keepParticle = false;
}


void Foam::trackedParticle::hitWallPatch
(
    Cloud<trackedParticle>&,
    trackingData& td
)
{
    // Remove particle
    td.keepParticle = false;
}


void Foam::trackedParticle::correctAfterParallelTransfer
(
    Cloud<trackedParticle>& cloud,
    trackingData& td
)
{
    particle::correctAfterParallelTransfer(cloud, td);

    // Mark edge we are currently on (if any). This was set on the sending
    // processor but has not yet been set on the receiving side.
    const label feati = i(), edgei = k();
    if (edgei != -1)
    {
        td.featureEdgeVisited_[feati].set(edgei, 1u);
    }
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<(Ostream& os, const trackedParticle& p)
{
    if (os.format() == IOstream::ASCII)
    {
        os  << static_cast<const particle&>(p)
            << token::SPACE << p.start_
            << token::SPACE << p.end_
            << token::SPACE << p.level_
            << token::SPACE << p.i_
            << token::SPACE << p.j_
            << token::SPACE << p.k_;
    }
    else
    {
        os  << static_cast<const particle&>(p);
        os.write
        (
            reinterpret_cast<const char*>(&p.start_),
            sizeof(p.start_) + sizeof(p.end_) + sizeof(p.level_)
          + sizeof(p.i_) + sizeof(p.j_) + sizeof(p.k_)
        );
    }

    // Check state of Ostream
    os.check("Ostream& operator<<(Ostream&, const trackedParticle&)");

    return os;
}


// ************************************************************************* //
