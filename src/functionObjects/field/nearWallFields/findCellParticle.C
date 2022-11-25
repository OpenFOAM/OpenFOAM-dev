/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2013-2022 OpenFOAM Foundation
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
    const label celli,
    const vector& displacement,
    const label data
)
:
    particle(mesh, position, celli),
    displacement_(displacement),
    data_(data)
{}


Foam::findCellParticle::findCellParticle(Istream& is, bool readFields)
:
    particle(is, readFields)
{
    if (readFields)
    {
        if (is.format() == IOstream::ASCII)
        {
            is >> displacement_;
            data_ = readLabel(is);
        }
        else
        {
            is.read
            (
                reinterpret_cast<char*>(&displacement_),
                sizeof(displacement_) + sizeof(data_)
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
    Cloud<findCellParticle>& cloud,
    trackingData& td
)
{
    td.keepParticle = true;
    td.sendToProc = -1;

    while (td.keepParticle && td.sendToProc == -1 && stepFraction() < 1)
    {
        const scalar f = 1 - stepFraction();
        trackToAndHitFace(f*displacement_, f, cloud, td);
    }

    if (td.sendToProc == -1)
    {
        // Hit endpoint or patch. If patch hit could do fancy stuff but just
        // to use the patch point is good enough for now.
        td.cellToData()[cell()].append(data());
        td.cellToEnd()[cell()].append(position(td.mesh));
    }

    return td.keepParticle;
}


void Foam::findCellParticle::hitWedgePatch
(
    Cloud<findCellParticle>&,
    trackingData& td
)
{
    // Remove particle
    td.keepParticle = false;
}


void Foam::findCellParticle::hitSymmetryPlanePatch
(
    Cloud<findCellParticle>&,
    trackingData& td
)
{
    // Remove particle
    td.keepParticle = false;
}


void Foam::findCellParticle::hitSymmetryPatch
(
    Cloud<findCellParticle>&,
    trackingData& td
)
{
    // Remove particle
    td.keepParticle = false;
}


void Foam::findCellParticle::hitCyclicPatch
(
    Cloud<findCellParticle>&,
    trackingData& td
)
{
    // Remove particle
    td.keepParticle = false;
}


void Foam::findCellParticle::hitProcessorPatch
(
    Cloud<findCellParticle>& cloud,
    trackingData& td
)
{
    const processorPolyPatch& ppp =
        static_cast<const processorPolyPatch&>
        (
            td.mesh.boundaryMesh()[patch(td.mesh)]
        );

    if (ppp.transform().transforms())
    {
        td.keepParticle = false;
    }
    else
    {
        particle::hitProcessorPatch(cloud, td);
    }
}


void Foam::findCellParticle::hitWallPatch
(
    Cloud<findCellParticle>&,
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
            << token::SPACE << p.displacement_
            << token::SPACE << p.data_;
    }
    else
    {
        os  << static_cast<const particle&>(p);
        os.write
        (
            reinterpret_cast<const char*>(&p.displacement_),
            sizeof(p.displacement_) + sizeof(p.data_)
        );
    }

    // Check state of Ostream
    os.check("Ostream& operator<<(Ostream&, const findCellParticle&)");

    return os;
}


// ************************************************************************* //
