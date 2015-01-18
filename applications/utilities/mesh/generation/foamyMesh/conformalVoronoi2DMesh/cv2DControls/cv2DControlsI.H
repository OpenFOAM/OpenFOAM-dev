/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013 OpenFOAM Foundation
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

inline Foam::scalar Foam::cv2DControls::minCellSize() const
{
    return minCellSize_;
}


inline Foam::scalar Foam::cv2DControls::minCellSize2() const
{
    return minCellSize2_;
}


inline Foam::scalar Foam::cv2DControls::maxQuadAngle() const
{
    return maxQuadAngle_;
}


inline Foam::scalar Foam::cv2DControls::nearWallAlignedDist() const
{
    return nearWallAlignedDist_;
}


inline Foam::scalar Foam::cv2DControls::nearWallAlignedDist2() const
{
    return nearWallAlignedDist2_;
}


inline Foam::Switch Foam::cv2DControls::insertSurfaceNearestPointPairs() const
{
    return insertSurfaceNearestPointPairs_;
}


inline Foam::Switch Foam::cv2DControls::mirrorPoints() const
{
    return mirrorPoints_;
}


inline Foam::Switch Foam::cv2DControls::insertSurfaceNearPointPairs() const
{
    return insertSurfaceNearPointPairs_;
}


inline Foam::Switch Foam::cv2DControls::objOutput() const
{
    return objOutput_;
}


inline Foam::Switch Foam::cv2DControls::meshedSurfaceOutput() const
{
    return meshedSurfaceOutput_;
}


inline Foam::Switch Foam::cv2DControls::randomiseInitialGrid() const
{
    return randomiseInitialGrid_;
}


inline Foam::scalar Foam::cv2DControls::randomPerturbation() const
{
    return randomPerturbation_;
}


inline Foam::label Foam::cv2DControls::maxBoundaryConformingIter() const
{
    return maxBoundaryConformingIter_;
}


inline Foam::scalar Foam::cv2DControls::span() const
{
    return span_;
}


inline Foam::scalar Foam::cv2DControls::span2() const
{
    return span2_;
}


inline Foam::scalar Foam::cv2DControls::minEdgeLen() const
{
    return minEdgeLen_;
}


inline Foam::scalar Foam::cv2DControls::minEdgeLen2() const
{
    return minEdgeLen2_;
}


inline Foam::scalar Foam::cv2DControls::maxNotchLen() const
{
    return maxNotchLen_;
}


inline Foam::scalar Foam::cv2DControls::maxNotchLen2() const
{
    return maxNotchLen2_;
}


inline Foam::scalar Foam::cv2DControls::minNearPointDist() const
{
    return minNearPointDist_;
}


inline Foam::scalar Foam::cv2DControls::minNearPointDist2() const
{
    return minNearPointDist2_;
}


inline Foam::scalar Foam::cv2DControls::ppDist() const
{
    return ppDist_;
}


// ************************************************************************* //
