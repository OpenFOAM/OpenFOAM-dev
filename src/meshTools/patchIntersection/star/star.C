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

#include "star.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::star::reset()
{
    forAll(starFaceFaces_, starFacei)
    {
        const label facei = starFaceFaces_[starFacei];
        if (facei != -1)
        {
            faceStarFaces_[starFaceFaces_[starFacei]] = -1;
        }
    }
    starFaceFaces_.clear();

    forAll(starEdgeEdges_, starEdgei)
    {
        const label edgei = starEdgeEdges_[starEdgei].edgei_;
        if (edgei != -1)
        {
            edgeStarEdges_[edgei] = -1;
        }
    }
    starEdgeEdges_.clear();
}


void Foam::star::swapStarEdges(const label starEdgeiA, const label starEdgeiB)
{
    const label starEdgeiA0 = starEdgeEdges_[starEdgeiA].starEdgei0_;
    const label starEdgeiA1 = starEdgeEdges_[starEdgeiA].starEdgei1_;

    const label starEdgeiB0 = starEdgeEdges_[starEdgeiB].starEdgei0_;
    const label starEdgeiB1 = starEdgeEdges_[starEdgeiB].starEdgei1_;

    if (starEdgeiA0 != -1)
    {
        starEdgeEdges_[starEdgeiA0].starEdgei1_ = starEdgeiB;
    }
    if (starEdgeiA1 != -1)
    {
        starEdgeEdges_[starEdgeiA1].starEdgei0_ = starEdgeiB;
    }

    if (starEdgeiB0 != -1)
    {
        starEdgeEdges_[starEdgeiB0].starEdgei1_ = starEdgeiA;
    }
    if (starEdgeiB1 != -1)
    {
        starEdgeEdges_[starEdgeiB1].starEdgei0_ = starEdgeiA;
    }

    Swap(starEdgeEdges_[starEdgeiA], starEdgeEdges_[starEdgeiB]);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::star::star()
:
    starFaceFaces_(),
    faceStarFaces_(),
    starEdgeEdges_(),
    edgeStarEdges_(),
    work_()
{}


Foam::star::context::context(star& s)
:
    star_(s)
{}


// * * * * * * * * * * * * * * * * Destructors * * * * * * * * * * * * * * * //

Foam::star::~star()
{}


Foam::star::context::~context()
{
    star_.reset();
}


// ************************************************************************* //
