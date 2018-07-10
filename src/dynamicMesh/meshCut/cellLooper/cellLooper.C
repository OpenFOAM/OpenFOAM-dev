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

#include "cellLooper.H"
#include "polyMesh.H"
#include "ListOps.H"
#include "meshTools.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(cellLooper, 0);
    defineRunTimeSelectionTable(cellLooper, word);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::cellLooper> Foam::cellLooper::New
(
    const word& type,
    const polyMesh& mesh
)
{
    wordConstructorTable::iterator cstrIter =
        wordConstructorTablePtr_->find(type);

    if (cstrIter == wordConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown set type "
            << type << nl << nl
            << "Valid cellLooper types : " << endl
            << wordConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<cellLooper>(cstrIter()(mesh));
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

Foam::labelList Foam::cellLooper::getVertFacesNonEdge
(
    const label celli,
    const label edgeI,
    const label vertI
) const
{
    // Get faces connected to startEdge
    label face0, face1;
    meshTools::getEdgeFaces(mesh(), celli, edgeI, face0, face1);

    const labelList& pFaces = mesh().pointFaces()[vertI];

    labelList vertFaces(pFaces.size());
    label vertFacei = 0;

    forAll(pFaces, pFacei)
    {
        label facei = pFaces[pFacei];

        if
        (
            (facei != face0)
         && (facei != face1)
         && (meshTools::faceOnCell(mesh(), celli, facei))
        )
        {
            vertFaces[vertFacei++] = facei;
        }
    }
    vertFaces.setSize(vertFacei);

    return vertFaces;
}


Foam::label Foam::cellLooper::getFirstVertEdge
(
    const label facei,
    const label vertI
) const
{
    const labelList& fEdges = mesh().faceEdges()[facei];

    forAll(fEdges, fEdgeI)
    {
        label edgeI = fEdges[fEdgeI];

        const edge& e = mesh().edges()[edgeI];

        if ((e.start() == vertI) || (e.end() == vertI))
        {
            return edgeI;
        }
    }

    FatalErrorInFunction
        << "Can not find edge on face " << facei
        << " using vertex " << vertI
        << abort(FatalError);

    return -1;
}


Foam::labelList Foam::cellLooper::getVertEdgesNonFace
(
    const label celli,
    const label facei,
    const label vertI
) const
{
    const labelList& exclEdges = mesh().faceEdges()[facei];

    const labelList& pEdges = mesh().pointEdges()[vertI];

    labelList vertEdges(pEdges.size());
    label vertEdgeI = 0;

    forAll(pEdges, pEdgeI)
    {
        label edgeI = pEdges[pEdgeI];

        if
        (
            (findIndex(exclEdges, edgeI) == -1)
         && meshTools::edgeOnCell(mesh(), celli, edgeI)
        )
        {
            vertEdges[vertEdgeI++] = edgeI;
        }
    }

    vertEdges.setSize(vertEdgeI);

    return vertEdges;
}


Foam::label Foam::cellLooper::getMisAlignedEdge
(
    const vector& refDir,
    const label celli
) const
{
    const labelList& cEdges = mesh().cellEdges()[celli];

    label cutEdgeI = -1;
    scalar maxCos = -great;

    forAll(cEdges, cEdgeI)
    {
        label edgeI = cEdges[cEdgeI];

        scalar cosAngle = mag(refDir & meshTools::normEdgeVec(mesh(), edgeI));

        if (cosAngle > maxCos)
        {
            maxCos = cosAngle;

            cutEdgeI = edgeI;
        }
    }

    return cutEdgeI;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::cellLooper::cellLooper(const polyMesh& mesh)
:
    edgeVertex(mesh)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::cellLooper::~cellLooper()
{}


// ************************************************************************* //
