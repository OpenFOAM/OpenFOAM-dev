/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
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

#include "ignitionSite.H"
#include "Time.H"
#include "volFields.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::ignitionSite::findIgnitionCells(const fvMesh& mesh)
{
    // Bit tricky: generate C and V before shortcutting if cannot find
    // cell locally. mesh.C generation uses parallel communication.
    const volVectorField& centres = mesh.C();
    const scalarField& vols = mesh.V();

    label ignCell = mesh.findCell(location_);
    if (ignCell == -1)
    {
        return;
    }

    scalar radius = diameter_/2.0;

    cells_.setSize(1);
    cellVolumes_.setSize(1);

    cells_[0] = ignCell;
    cellVolumes_[0] = vols[ignCell];

    scalar minDist = great;
    label nIgnCells = 1;

    forAll(centres, celli)
    {
        scalar dist = mag(centres[celli] - location_);

        if (dist < minDist)
        {
            minDist = dist;
        }

        if (dist < radius && celli != ignCell)
        {
            cells_.setSize(nIgnCells+1);
            cellVolumes_.setSize(nIgnCells+1);

            cells_[nIgnCells] = celli;
            cellVolumes_[nIgnCells] = vols[celli];

            nIgnCells++;
        }
    }

    if (cells_.size())
    {
        Pout<< "Found ignition cells:" << endl << cells_ << endl;
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::labelList& Foam::ignitionSite::cells() const
{
    if (mesh_.changing() && timeIndex_ != db_.timeIndex())
    {
        const_cast<ignitionSite&>(*this).findIgnitionCells(mesh_);
    }
    timeIndex_ = db_.timeIndex();

    return cells_;
}


bool Foam::ignitionSite::igniting() const
{
    scalar curTime = db_.value();
    scalar deltaT = db_.deltaTValue();

    return
    (
        (curTime - deltaT >= time_)
        &&
        (curTime - deltaT < time_ + max(duration_, deltaT) + small)
    );
}


bool Foam::ignitionSite::ignited() const
{
    scalar curTime = db_.value();
    scalar deltaT = db_.deltaTValue();

    return(curTime - deltaT >= time_);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

void Foam::ignitionSite::operator=(const ignitionSite& is)
{
    location_ = is.location_;
    diameter_ = is.diameter_;
    time_ = is.time_;
    duration_ = is.duration_;
    strength_ = is.strength_;
    cells_ = is.cells_;
    cellVolumes_ = is.cellVolumes_;
}


// ************************************************************************* //
