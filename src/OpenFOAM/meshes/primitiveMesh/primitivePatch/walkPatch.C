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

#include "walkPatch.H"
#include "ListOps.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(walkPatch, 0);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::label Foam::walkPatch::getNeighbour
(
    const label facei,
    const label fp,
    const label v0,
    const label v1
) const
{
    const labelList& fEdges = pp_.faceEdges()[facei];

    const edgeList& edges = pp_.edges();


    label nbrEdgeI = -1;

    // Shortcut: maybe faceEdges are sorted(?) in which case fEdges[fp] is
    // edge between v0 and v1.
    const edge& e = edges[fEdges[fp]];

    if ((e[0] == v0 && e[1] == v1) || (e[0] == v1 && e[1] == v0))
    {
        // Correct edge.
        nbrEdgeI = fEdges[fp];
    }
    else
    {
        // Loop over all faceEdges.
        forAll(fEdges, i)
        {
            label edgeI = fEdges[i];

            const edge& e = edges[edgeI];

            if
            (
                (e[0] == v0 && e[1] == v1)
             || (e[0] == v1 && e[1] == v0)
            )
            {
                // Found edge on face which uses v0, v1.
                nbrEdgeI = edgeI;

                break;
            }
        }
    }


    if (nbrEdgeI == -1)
    {
        FatalErrorInFunction
            << "Did not find edge on face " << facei << " that uses vertices"
            << v0 << " and " << v1 << abort(FatalError);
    }


    // Get neighbouring face.

    const labelList& eFaces = pp_.edgeFaces()[nbrEdgeI];

    if (eFaces.size() == 1)
    {
        return -1;
    }
    else if (eFaces.size() == 2)
    {
        label nbrFacei = eFaces[0];

        if (nbrFacei == facei)
        {
            nbrFacei = eFaces[1];
        }

        return nbrFacei;
    }
    else
    {
        FatalErrorInFunction
            << "Illegal surface on patch. Face " << facei
            << " at vertices " << v0 << ',' << v1
            << " has fewer than 1 or more than 2 neighbours"
            << abort(FatalError);
        return -1;
    }
}


void Foam::walkPatch::faceToFace
(
    const labelList& changedFaces,
    const labelList& enterVerts,

    labelList& nbrFaces,
    labelList& nbrEnterVerts
)
{
    nbrFaces.setSize(pp_.size());
    nbrEnterVerts.setSize(pp_.size());
    label changedI = 0;

    forAll(changedFaces, i)
    {
        label facei = changedFaces[i];
        label enterVertI = enterVerts[i];

        if (!visited_[facei])
        {
            // Do this face
            visited_[facei] = true;
            visitOrder_.append(facei);

            const face& f = pp_.localFaces()[facei];

            label fp = findIndex(f, enterVertI);

            indexInFace_.append(fp);

            // Visit neighbouring faces in order, starting at fp.
            forAll(f, i)
            {
                label fp1 = reverse_ ? f.rcIndex(fp) : f.fcIndex(fp);
                label nbr = getNeighbour(facei, fp, f[fp], f[fp1]);

                if
                (
                    nbr != -1
                 && !visited_[nbr]
                 && faceZone_[nbr] == faceZone_[facei]
                )
                {
                    nbrFaces[changedI] = nbr;
                    nbrEnterVerts[changedI] = f[fp];
                    changedI++;
                }

                fp = fp1;
            }
        }
    }

    nbrFaces.setSize(changedI);
    nbrEnterVerts.setSize(changedI);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::walkPatch::walkPatch
(
    const primitivePatch& pp,
    const labelList& faceZone,
    const bool reverse,
    const label facei,
    const label enterVertI,
    boolList& visited
)
:
    pp_(pp),
    faceZone_(faceZone),
    reverse_(reverse),
    visited_(visited),
    visitOrder_(pp.size()),
    indexInFace_(pp.size())
{
    // List of faces that have been visited in the current iteration.
    labelList changedFaces(1, facei);
    // Corresponding list of entry vertices
    labelList enterVerts(1, enterVertI);

    while (true)
    {
        labelList nbrFaces;
        labelList nbrEnterVerts;

        faceToFace
        (
            changedFaces,
            enterVerts,

            nbrFaces,
            nbrEnterVerts
        );


        if (nbrFaces.empty())
        {
            break;
        }

        changedFaces = nbrFaces;
        enterVerts = nbrEnterVerts;
    }

    visitOrder_.shrink();
    indexInFace_.shrink();
}


// ************************************************************************* //
