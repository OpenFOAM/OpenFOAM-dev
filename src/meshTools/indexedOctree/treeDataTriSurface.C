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

#include "treeDataTriSurface.H"
#include "triSurfaceTools.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<>
Foam::volumeType Foam::treeDataPrimitivePatch<Foam::triSurface>::getVolumeType
(
    const indexedOctree<treeDataPrimitivePatch<triSurface>>& oc,
    const point& sample
) const
{
    // Find nearest face to sample
    pointIndexHit info = oc.findNearest(sample, sqr(great));

    if (info.index() == -1)
    {
        FatalErrorInFunction
            << "Could not find " << sample << " in octree."
            << abort(FatalError);
    }

    // Get actual intersection point on face
    label facei = info.index();

    triSurfaceTools::sideType t = triSurfaceTools::surfaceSide
    (
        patch_,
        sample,
        facei
    );

    if (t == triSurfaceTools::UNKNOWN)
    {
        return volumeType::unknown;
    }
    else if (t == triSurfaceTools::INSIDE)
    {
        return volumeType::inside;
    }
    else if (t == triSurfaceTools::OUTSIDE)
    {
        return volumeType::outside;
    }
    else
    {
        FatalErrorInFunction
            << "problem" << abort(FatalError);
        return volumeType::unknown;
    }
}


// ************************************************************************* //
