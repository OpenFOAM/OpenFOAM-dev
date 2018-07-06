/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2012-2018 OpenFOAM Foundation
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

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::searchableSurfaces& Foam::conformationSurfaces::geometry() const
{
    return allGeometry_;
}


const Foam::PtrList<Foam::extendedFeatureEdgeMesh>&
Foam::conformationSurfaces::features() const
{
    return features_;
}


const Foam::point& Foam::conformationSurfaces::locationInMesh() const
{
    return locationInMesh_;
}


const Foam::labelList& Foam::conformationSurfaces::surfaces() const
{
    return surfaces_;
}


const Foam::List<Foam::word>& Foam::conformationSurfaces::patchNames() const
{
    return patchNames_;
}


const Foam::PtrList<Foam::surfaceZonesInfo>&
Foam::conformationSurfaces::surfZones() const
{
    return surfZones_;
}


const Foam::PtrList<Foam::dictionary>&
Foam::conformationSurfaces::patchInfo() const
{
    return patchInfo_;
}


const Foam::treeBoundBox& Foam::conformationSurfaces::globalBounds() const
{
    return globalBounds_;
}


// ************************************************************************* //
