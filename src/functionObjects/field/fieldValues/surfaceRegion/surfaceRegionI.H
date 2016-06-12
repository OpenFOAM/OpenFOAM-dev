/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
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

#include "surfaceRegion.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline const Foam::functionObjects::fieldValues::surfaceRegion::regionTypes&
Foam::functionObjects::fieldValues::surfaceRegion::regionType() const
{
    return regionType_;
}


inline const Foam::labelList&
Foam::functionObjects::fieldValues::surfaceRegion::faceId() const
{
    return faceId_;
}


inline const Foam::labelList&
Foam::functionObjects::fieldValues::surfaceRegion::facePatch() const
{
    return facePatchId_;
}


inline const Foam::labelList&
Foam::functionObjects::fieldValues::surfaceRegion::faceSign() const
{
    return faceSign_;
}


// ************************************************************************* //
