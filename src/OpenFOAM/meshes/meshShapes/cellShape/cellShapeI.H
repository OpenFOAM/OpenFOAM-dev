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

#include "Istream.H"
#include "cell.H"
#include "cellModeller.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline Foam::cellShape::cellShape()
:
    m(nullptr)
{}


inline Foam::cellShape::cellShape
(
    const cellModel& M,
    const labelList& l,
    const bool doCollapse
)
:
    labelList(l),
    m(&M)
{
    if (doCollapse)
    {
        collapse();
    }
}


inline Foam::cellShape::cellShape
(
    const word& model,
    const labelList& l,
    const bool doCollapse
)
:
    labelList(l),
    m(cellModeller::lookup(model))
{
    if (doCollapse)
    {
        collapse();
    }
}


inline Foam::cellShape::cellShape(Istream& is)
{
    is >> *this;
}


inline Foam::autoPtr<Foam::cellShape> Foam::cellShape::clone() const
{
    return autoPtr<cellShape>(new cellShape(*this));
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline Foam::pointField Foam::cellShape::points
(
    const pointField& meshPoints
) const
{
    // There are as many points as there labels for them
    pointField p(size());

    // For each point in list, set it to the point in 'pnts' addressed
    // by 'labs'
    forAll(p, i)
    {
        p[i] = meshPoints[operator[](i)];
    }

    // Return list
    return p;
}


inline const Foam::cellModel& Foam::cellShape::model() const
{
    return *m;
}


inline Foam::labelList Foam::cellShape::meshFaces
(
    const faceList& allFaces,
    const cell& cFaces
) const
{
    // Faces in model order
    faceList localFaces(faces());

    // Do linear match (usually cell shape is low complexity)

    labelList modelToMesh(localFaces.size(), -1);

    forAll(localFaces, i)
    {
        const face& localF = localFaces[i];

        forAll(cFaces, j)
        {
            label meshFacei = cFaces[j];

            if (allFaces[meshFacei] == localF)
            {
                modelToMesh[i] = meshFacei;

                break;
            }
        }
    }

    return modelToMesh;
}


inline Foam::labelList Foam::cellShape::meshEdges
(
    const edgeList& allEdges,
    const labelList& cEdges
) const
{
    // Edges in model order
    edgeList localEdges(edges());

    // Do linear match (usually cell shape is low complexity)

    labelList modelToMesh(localEdges.size(), -1);

    forAll(localEdges, i)
    {
        const edge& e = localEdges[i];

        forAll(cEdges, j)
        {
            label edgeI = cEdges[j];

            if (allEdges[edgeI] == e)
            {
                modelToMesh[i] = edgeI;

                break;
            }
        }
    }

    return modelToMesh;
}


inline Foam::faceList Foam::cellShape::faces() const
{
    return m->faces(*this);
}


inline Foam::faceList Foam::cellShape::collapsedFaces() const
{
    faceList oldFaces(faces());

    faceList newFaces(oldFaces.size());
    label newFacei = 0;

    forAll(oldFaces, oldFacei)
    {
        const face& f = oldFaces[oldFacei];

        face& newF = newFaces[newFacei];

        newF.setSize(f.size());

        label newFp = 0;
        label prevVertI = -1;

        forAll(f, fp)
        {
            label vertI = f[fp];

            if (vertI != prevVertI)
            {
                newF[newFp++] = vertI;

                prevVertI = vertI;
            }
        }

        if ((newFp > 1) && (newF[newFp-1] == newF[0]))
        {
            --newFp;
        }

        if (newFp > 2)
        {
            // Size face and go to next one
            newF.setSize(newFp);

            newFacei++;
        }
    }
    newFaces.setSize(newFacei);

    return newFaces;
}


inline Foam::label Foam::cellShape::nFaces() const
{
    return m->nFaces();
}


inline Foam::edgeList Foam::cellShape::edges() const
{
    return m->edges(*this);
}


inline Foam::label Foam::cellShape::nEdges() const
{
    return m->nEdges();
}


inline Foam::label Foam::cellShape::nPoints() const
{
    return size();
}


inline Foam::point Foam::cellShape::centre(const pointField& points) const
{
    return m->centre(*this, points);
}


inline Foam::scalar Foam::cellShape::mag(const pointField& points) const
{
    return m->mag(*this, points);
}


// ************************************************************************* //
