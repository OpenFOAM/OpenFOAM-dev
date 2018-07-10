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

#include "patchZones.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(patchZones, 0);
}


// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

Foam::labelList Foam::patchZones::faceToEdge
(
    const labelList& changedFaces,
    labelList& edgeRegion
)
{
    labelList changedEdges(pp_.nEdges(), -1);
    label changedI = 0;

    forAll(changedFaces, i)
    {
        label facei = changedFaces[i];

        const labelList& fEdges = pp_.faceEdges()[facei];

        forAll(fEdges, fEdgeI)
        {
            label edgeI = fEdges[fEdgeI];

            if (!borderEdge_[edgeI] && (edgeRegion[edgeI] == -1))
            {
                edgeRegion[edgeI] = nZones_;

                changedEdges[changedI++] = edgeI;
            }
        }
    }

    changedEdges.setSize(changedI);

    return changedEdges;
}


Foam::labelList Foam::patchZones::edgeToFace(const labelList& changedEdges)
{
    labelList changedFaces(pp_.size(), -1);
    label changedI = 0;

    forAll(changedEdges, i)
    {
        label edgeI = changedEdges[i];

        const labelList& eFaces = pp_.edgeFaces()[edgeI];

        forAll(eFaces, eFacei)
        {
            label facei = eFaces[eFacei];

            if (operator[](facei) == -1)
            {
                operator[](facei) = nZones_;

                changedFaces[changedI++] = facei;
            }
        }
    }

    changedFaces.setSize(changedI);

    return changedFaces;
}


void Foam::patchZones::markZone(label facei)
{
    // List of faces whose faceZone has been set.
    labelList changedFaces(1, facei);
    // List of edges whose faceZone has been set.
    labelList changedEdges;

    // Zones on all edges.
    labelList edgeZone(pp_.nEdges(), -1);

    while (true)
    {
        changedEdges = faceToEdge(changedFaces, edgeZone);

        if (debug)
        {
            Info<< "From changedFaces:" << changedFaces.size()
                << " to changedEdges:" << changedEdges.size()
                << endl;
        }

        if (changedEdges.empty())
        {
            break;
        }

        changedFaces = edgeToFace(changedEdges);

        if (debug)
        {
            Info<< "From changedEdges:" << changedEdges.size()
                << " to changedFaces:" << changedFaces.size()
                << endl;
        }

        if (changedEdges.empty())
        {
            break;
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::patchZones::patchZones
(
    const primitivePatch& pp,
    const boolList& borderEdge
)
:
    labelList(pp.size(), -1),
    pp_(pp),
    borderEdge_(borderEdge),
    nZones_(0)
{
    // Finds areas delimited by borderEdge (or 'real' edges).
    // Fills *this with zone number accordingly.

    if (borderEdge.size() != pp_.nEdges())
    {
        FatalErrorInFunction
            << "borderEdge boolList not same size as number of edges" << endl
            << "borderEdge:" << borderEdge.size() << endl
            << "nEdges    :" << pp_.nEdges()
            << abort(FatalError);
    }

    label facei = 0;

    while (true)
    {
        // Find first non-visited face
        for (; facei < pp_.size(); facei++)
        {
            if (operator[](facei) == -1)
            {
                operator[](facei) = nZones_;

                markZone(facei);

                break;
            }
        }

        if (facei == pp_.size())
        {
            // Finished.
            break;
        }

        nZones_++;
    }
}


// ************************************************************************* //
