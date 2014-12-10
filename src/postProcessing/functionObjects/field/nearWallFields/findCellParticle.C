/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013 OpenFOAM Foundation
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

#include "findCellParticle.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::findCellParticle::findCellParticle
(
    const polyMesh& mesh,
    const vector& position,
    const label cellI,
    const label tetFaceI,
    const label tetPtI,
    const point& end,
    const label data
)
:
    particle(mesh, position, cellI, tetFaceI, tetPtI),
    end_(end),
    data_(data)
{}


Foam::findCellParticle::findCellParticle
(
    const polyMesh& mesh,
    Istream& is,
    bool readFields
)
:
    particle(mesh, is, readFields)
{
    if (readFields)
    {
        if (is.format() == IOstream::ASCII)
        {
            is >> end_;
            data_ = readLabel(is);
        }
        else
        {
            is.read
            (
                reinterpret_cast<char*>(&end_),
                sizeof(end_) + sizeof(data_)
            );
        }
    }

    // Check state of Istream
    is.check
    (
        "findCellParticle::findCellParticle"
        "(const Cloud<findCellParticle>&, Istream&, bool)"
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::findCellParticle::move
(
    trackingData& td,
    const scalar maxTrackLen
)
{
    td.switchProcessor = false;
    td.keepParticle = true;

    scalar tEnd = (1.0 - stepFraction())*maxTrackLen;
    scalar dtMax = tEnd;

    while (td.keepParticle && !td.switchProcessor && tEnd > SMALL)
    {
        // set the lagrangian time-step
        scalar dt = min(dtMax, tEnd);

        dt *= trackToFace(end_, td);

        tEnd -= dt;
        stepFraction() = 1.0 - tEnd/maxTrackLen;
    }

    if (tEnd < SMALL || !td.keepParticle)
    {
        // Hit endpoint or patch. If patch hit could do fancy stuff but just
        // to use the patch point is good enough for now.
        td.cellToData()[cell()].append(data());
        td.cellToEnd()[cell()].append(position());
    }

    return td.keepParticle;
}


bool Foam::findCellParticle::hitPatch
(
    const polyPatch&,
    trackingData& td,
    const label patchI,
    const scalar trackFraction,
    const tetIndices& tetIs
)
{
    return false;
}


void Foam::findCellParticle::hitWedgePatch
(
    const wedgePolyPatch&,
    trackingData& td
)
{
    // Remove particle
    td.keepParticle = false;
}


void Foam::findCellParticle::hitSymmetryPlanePatch
(
    const symmetryPlanePolyPatch&,
    trackingData& td
)
{
    // Remove particle
    td.keepParticle = false;
}


void Foam::findCellParticle::hitSymmetryPatch
(
    const symmetryPolyPatch&,
    trackingData& td
)
{
    // Remove particle
    td.keepParticle = false;
}


void Foam::findCellParticle::hitCyclicPatch
(
    const cyclicPolyPatch&,
    trackingData& td
)
{
    // Remove particle
    td.keepParticle = false;
}


void Foam::findCellParticle::hitProcessorPatch
(
    const processorPolyPatch&,
    trackingData& td
)
{
    // Remove particle
    td.switchProcessor = true;
}


void Foam::findCellParticle::hitWallPatch
(
    const wallPolyPatch& wpp,
    trackingData& td,
    const tetIndices&
)
{
    // Remove particle
    td.keepParticle = false;
}


void Foam::findCellParticle::hitPatch
(
    const polyPatch& wpp,
    trackingData& td
)
{
    // Remove particle
    td.keepParticle = false;
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<(Ostream& os, const findCellParticle& p)
{
    if (os.format() == IOstream::ASCII)
    {
        os  << static_cast<const particle&>(p)
            << token::SPACE << p.end_
            << token::SPACE << p.data_;
    }
    else
    {
        os  << static_cast<const particle&>(p);
        os.write
        (
            reinterpret_cast<const char*>(&p.end_),
            sizeof(p.end_) + sizeof(p.data_)
        );
    }

    // Check state of Ostream
    os.check("Ostream& operator<<(Ostream&, const findCellParticle&)");

    return os;
}


// ************************************************************************* //
