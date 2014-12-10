/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

#include "PatchTools.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template
<
    class Face,
    template<class> class FaceList,
    class PointField,
    class PointType
>

Foam::labelList
Foam::PatchTools::edgeOwner
(
    const PrimitivePatch<Face, FaceList, PointField, PointType>& p
)
{
    const edgeList& edges = p.edges();
    const labelListList& edgeFaces = p.edgeFaces();
    const List<Face>& localFaces = p.localFaces();

    // create the owner list
    labelList edgeOwner(edges.size(), -1);

    forAll(edges, edgeI)
    {
        const labelList& nbrFaces = edgeFaces[edgeI];

        if (nbrFaces.size() == 1)
        {
            edgeOwner[edgeI] = nbrFaces[0];
        }
        else
        {
            // Find the first face whose vertices are aligned with the edge.
            // with multiply connected edges, this is the best we can do
            forAll(nbrFaces, i)
            {
                const Face& f = localFaces[nbrFaces[i]];

                if (f.edgeDirection(edges[edgeI]) > 0)
                {
                    edgeOwner[edgeI] = nbrFaces[i];
                    break;
                }
            }

            if (edgeOwner[edgeI] == -1)
            {
                FatalErrorIn
                (
                    "PatchTools::edgeOwner()"
                )
                    << "Edge " << edgeI << " vertices:" << edges[edgeI]
                    << " is used by faces " << nbrFaces
                    << " vertices:"
                    << UIndirectList<Face>(localFaces, nbrFaces)()
                    << " none of which use the edge vertices in the same order"
                    << nl << "I give up" << abort(FatalError);
            }
        }
    }

    return edgeOwner;
}


// ************************************************************************* //
