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

#include "polyMesh.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline const Foam::polyCellSet::selectionTypes&
Foam::polyCellSet::selectionType() const
{
    return selectionType_;
}


inline const Foam::word& Foam::polyCellSet::cellSetName() const
{
    return cellSetName_;
}


inline Foam::labelUList Foam::polyCellSet::cells() const
{
    if (selectionType_ == selectionTypes::all)
    {
        return identityMap(mesh_.nCells());
    }
    else if (selectionType_ == selectionTypes::cellZone)
    {
        return mesh_.cellZones()[cellSetName_];
    }
    else
    {
        return cells_;
    }
}


inline Foam::label Foam::polyCellSet::nCells() const
{
    if (selectionType_ == selectionTypes::all)
    {
        return mesh_.nCells();
    }
    else
    {
        return cells().size();
    }
}


inline bool Foam::polyCellSet::all() const
{
    return selectionType_ == selectionTypes::all;
}


inline Foam::label Foam::polyCellSet::celli(const label i) const
{
    if (selectionType_ == selectionTypes::all)
    {
        return i;
    }
    else if (selectionType_ == selectionTypes::cellZone)
    {
        return mesh_.cellZones()[cellSetName_][i];
    }
    else
    {
        return cells_[i];
    }
}


// ************************************************************************* //
