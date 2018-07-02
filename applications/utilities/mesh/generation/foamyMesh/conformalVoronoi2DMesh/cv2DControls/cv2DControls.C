/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013-2018 OpenFOAM Foundation
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

#include "cv2DControls.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


Foam::cv2DControls::cv2DControls
(
    const dictionary& controlDict,
    const boundBox& bb
)
:
    dict_(controlDict),

    motionControl_(controlDict.subDict("motionControl")),
    conformationControl_(controlDict.subDict("surfaceConformation")),

    minCellSize_(readScalar(motionControl_.lookup("minCellSize"))),
    minCellSize2_(Foam::sqr(minCellSize_)),

    maxQuadAngle_(readScalar(conformationControl_.lookup("maxQuadAngle"))),

    nearWallAlignedDist_
    (
        readScalar(motionControl_.lookup("nearWallAlignedDist"))*minCellSize_
    ),
    nearWallAlignedDist2_(Foam::sqr(nearWallAlignedDist_)),

    insertSurfaceNearestPointPairs_
    (
            conformationControl_.lookup("insertSurfaceNearestPointPairs")
    ),
    mirrorPoints_(conformationControl_.lookup("mirrorPoints")),
    insertSurfaceNearPointPairs_
    (
        conformationControl_.lookup("insertSurfaceNearPointPairs")
    ),

    objOutput_(motionControl_.lookupOrDefault<Switch>("objOutput", false)),

    meshedSurfaceOutput_
    (
        motionControl_.lookupOrDefault<Switch>("meshedSurfaceOutput", false)
    ),

    randomiseInitialGrid_(conformationControl_.lookup("randomiseInitialGrid")),
    randomPerturbation_
    (
        readScalar(conformationControl_.lookup("randomPerturbation"))
    ),

    maxBoundaryConformingIter_
    (
        readLabel(conformationControl_.lookup("maxBoundaryConformingIter"))
    ),

    span_
    (
        max(mag(bb.max().x()), mag(bb.min().x()))
      + max(mag(bb.max().y()), mag(bb.min().y()))
    ),
    span2_(Foam::sqr(span_)),

    minEdgeLen_
    (
        readScalar(conformationControl_.lookup("minEdgeLenCoeff"))
       *minCellSize_
    ),
    minEdgeLen2_(Foam::sqr(minEdgeLen_)),

    maxNotchLen_
    (
        readScalar(conformationControl_.lookup("maxNotchLenCoeff"))
       *minCellSize_
    ),
    maxNotchLen2_(Foam::sqr(maxNotchLen_)),

    minNearPointDist_
    (
        readScalar(conformationControl_.lookup("minNearPointDistCoeff"))
       *minCellSize_
    ),
    minNearPointDist2_(Foam::sqr(minNearPointDist_)),

    ppDist_
    (
        readScalar(conformationControl_.lookup("pointPairDistanceCoeff"))
       *minCellSize_
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::cv2DControls::~cv2DControls()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::cv2DControls::write(Ostream& os) const
{
    os.indentLevel() = 1;
    os.precision(2);
    os.flags(ios_base::scientific);

    os << nl << "Outputting CV2D Mesher controls:" << nl
       << token::BEGIN_BLOCK << nl
       << indent << "minCellSize2_         : " << minCellSize2_ << nl
       << indent << "span_ / span2_        : " << span_ << " / " << span2_ << nl
       << indent << "maxNotchLen2_         : " << maxNotchLen2_ << nl
       << indent << "minNearPointDist2_    : " << minNearPointDist2_ << nl
       << indent << "nearWallAlignedDist2_ : " << nearWallAlignedDist2_ << nl
       << indent << "ppDist_               : " << ppDist_ << nl
       << indent << "minEdgeLen2_          : " << minEdgeLen2_ << nl
       << token::END_BLOCK << endl;
}


// * * * * * * * * * * * * * * IOStream operators  * * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<(Ostream& os, const cv2DControls& s)
{
    s.write(os);
    return os;
}



// ************************************************************************* //
