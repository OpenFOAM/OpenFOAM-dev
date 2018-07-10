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

InClass
    vtkPVFoam

\*---------------------------------------------------------------------------*/

#ifndef vtkOpenFOAMPoints_H
#define vtkOpenFOAMPoints_H

// VTK includes
#include "vtkPoints.h"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

inline void vtkInsertNextOpenFOAMPoint
(
    vtkPoints *points,
    const Foam::point& p
)
{
    points->InsertNextPoint(p.x(), p.y(), p.z());
}

#if 0
// this should be faster, but didn't get it working ...
inline void vtkSetOpenFOAMPoint
(
    vtkPoints *points,
    const Foam::label id,
    const Foam::point& p
)
{
    points->SetPoint(id, p.x(), p.y(), p.z());
}


// Convert OpenFOAM mesh vertices to VTK
inline vtkPoints* vtkSetOpenFOAMPoints(const Foam::pointField& points)
{
    vtkPoints *vtkpoints = vtkPoints::New();
    vtkpoints->SetNumberOfPoints(points.size());
    forAll(points, i)
    {
        const Foam::point& p = points[i];
        vtkpoints->SetPoint(i, p.x(), p.y(), p.z());
    }

    return vtkpoints;
}

#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
