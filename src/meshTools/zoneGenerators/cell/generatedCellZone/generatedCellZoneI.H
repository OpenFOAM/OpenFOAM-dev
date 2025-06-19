/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2021-2025 OpenFOAM Foundation
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

#include "polyMesh.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline bool Foam::generatedCellZone::all() const
{
    return all_;
}


inline const Foam::word& Foam::generatedCellZone::zoneName() const
{
    static word allName_("all");

    if (all_)
    {
        return allName_;
    }
    else
    {
        return cellZone_.cZone().name();
    }
}


inline const Foam::cellZone& Foam::generatedCellZone::zone() const
{
    if (all_)
    {
        return mesh_.cellZones().all();
    }
    else
    {
        return cellZone_.cZone();
    }
}


inline Foam::label Foam::generatedCellZone::nCells() const
{
    if (all_)
    {
        return mesh_.nCells();
    }
    else
    {
        return cellZone_.cZone().size();
    }
}


inline Foam::label Foam::generatedCellZone::celli(const label i) const
{
    if (all_)
    {
        return i;
    }
    else
    {
        return cellZone_.cZone()[i];
    }
}


// ************************************************************************* //
