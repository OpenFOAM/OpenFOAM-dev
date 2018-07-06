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

#include "PrimitivePatch.H"
#include "HashSet.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template
<
    class Face,
    template<class> class FaceList,
    class PointField,
    class PointType
>
void
Foam::PrimitivePatch<Face, FaceList, PointField, PointType>::
calcBdryPoints() const
{
    if (debug)
    {
        InfoInFunction << "Calculating boundary points" << endl;
    }

    if (boundaryPointsPtr_)
    {
        // it is considered an error to attempt to recalculate
        // if already allocated
        FatalErrorInFunction
            << "edge types already calculated"
            << abort(FatalError);
    }

    const edgeList& e = edges();

    labelHashSet bp(2*e.size());

    for (label edgeI = nInternalEdges_; edgeI < e.size(); edgeI++)
    {
        const edge& curEdge = e[edgeI];

        bp.insert(curEdge.start());
        bp.insert(curEdge.end());
    }

    boundaryPointsPtr_ = new labelList(bp.toc());
    sort(*boundaryPointsPtr_);

    if (debug)
    {
        Info<< "    Finished." << endl;
    }
}


// ************************************************************************* //
