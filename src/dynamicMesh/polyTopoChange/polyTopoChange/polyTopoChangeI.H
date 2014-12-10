/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

#include "face.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

inline bool Foam::polyTopoChange::pointRemoved(const label pointI) const
{
    const point& pt = points_[pointI];

    return
        pt.x() > 0.5*vector::max.x()
     && pt.y() > 0.5*vector::max.y()
     && pt.z() > 0.5*vector::max.z();
}


inline bool Foam::polyTopoChange::faceRemoved(const label faceI) const
{
    return faces_[faceI].empty();
}


inline bool Foam::polyTopoChange::cellRemoved(const label cellI) const
{
    return cellMap_[cellI] == -2;
}


inline void Foam::polyTopoChange::setNumPatches(const label nPatches)
{
    nPatches_ = nPatches;
}


// ************************************************************************* //
