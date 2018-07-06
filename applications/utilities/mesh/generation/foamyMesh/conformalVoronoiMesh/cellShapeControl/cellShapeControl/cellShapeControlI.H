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

inline Foam::cellShapeControlMesh&
Foam::cellShapeControl::shapeControlMesh()
{
    return shapeControlMesh_;
}


inline const Foam::cellShapeControlMesh&
Foam::cellShapeControl::shapeControlMesh() const
{
    return shapeControlMesh_;
}


inline const Foam::scalar& Foam::cellShapeControl::defaultCellSize() const
{
    return defaultCellSize_;
}


inline const Foam::cellAspectRatioControl&
Foam::cellShapeControl::aspectRatio() const
{
    return aspectRatio_;
}


inline const Foam::cellSizeAndAlignmentControls&
Foam::cellShapeControl::sizeAndAlignment() const
{
    return sizeAndAlignment_;
}


inline const Foam::scalar& Foam::cellShapeControl::minimumCellSize() const
{
    return minimumCellSize_;
}


// ************************************************************************* //
