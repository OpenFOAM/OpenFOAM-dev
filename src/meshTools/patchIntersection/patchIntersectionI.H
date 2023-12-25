/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2023 OpenFOAM Foundation
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

#include "patchIntersection.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

inline const Foam::pointField& Foam::patchIntersection::points() const
{
    return points_;
}


inline const Foam::labelList&
Foam::patchIntersection::srcPointPoints() const
{
    return srcPointPoints_;
}


inline const Foam::labelList&
Foam::patchIntersection::tgtPointPoints() const
{
    return tgtPointPoints_;
}


inline const Foam::DynamicList<Foam::label>&
Foam::patchIntersection::pointSrcPoints() const
{
    return pointSrcPoints_;
}


inline const Foam::DynamicList<Foam::label>&
Foam::patchIntersection::pointTgtPoints() const
{
    return pointTgtPoints_;
}


inline const Foam::List<Foam::DynamicList<Foam::label>>&
Foam::patchIntersection::srcEdgePoints() const
{
    return srcEdgePoints_;
}


inline const Foam::List<Foam::DynamicList<Foam::label>>&
Foam::patchIntersection::tgtEdgePoints() const
{
    return tgtEdgePoints_;
}


inline const Foam::DynamicList<Foam::label>&
Foam::patchIntersection::pointSrcEdges() const
{
    return pointSrcEdges_;
}


inline const Foam::DynamicList<Foam::label>&
Foam::patchIntersection::pointTgtEdges() const
{
    return pointTgtEdges_;
}


inline const Foam::DynamicList<Foam::label>&
Foam::patchIntersection::pointSrcFaces() const
{
    return pointSrcFaces_;
}


inline const Foam::DynamicList<Foam::label>&
Foam::patchIntersection::pointTgtFaces() const
{
    return pointTgtFaces_;
}


inline const Foam::faceList& Foam::patchIntersection::faces() const
{
    return faces_;
}


inline const Foam::List<Foam::DynamicList<Foam::label>>&
Foam::patchIntersection::srcFaceFaces() const
{
    return srcFaceFaces_;
}


inline const Foam::List<Foam::DynamicList<Foam::label>>&
Foam::patchIntersection::tgtFaceFaces() const
{
    return tgtFaceFaces_;
}


inline const Foam::DynamicList<Foam::label>&
Foam::patchIntersection::faceSrcFaces() const
{
    return faceSrcFaces_;
}


inline const Foam::DynamicList<Foam::label>&
Foam::patchIntersection::faceTgtFaces() const
{
    return faceTgtFaces_;
}


// ************************************************************************* //
