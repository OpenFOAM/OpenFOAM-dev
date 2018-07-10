/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
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

#include "DelaunayMeshTools.H"
#include "meshTools.H"
#include "OFstream.H"
#include "pointConversion.H"
#include "pointIOField.H"
#include "indexedVertexOps.H"

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class Triangulation>
void Foam::DelaunayMeshTools::writeOBJ
(
    const fileName& fName,
    const Triangulation& t,
    const indexedVertexEnum::vertexType startPointType,
    const indexedVertexEnum::vertexType endPointType
)
{
    OFstream str(fName);

    Pout<< nl
        << "Writing points of types:" << nl;

    forAllConstIter
    (
        HashTable<unsigned int>,
        indexedVertexEnum::vertexTypeNames_,
        iter
    )
    {
        if (iter() >= startPointType && iter() <= endPointType)
        {
            Pout<< "    " << iter.key() << nl;
        }
    }

    Pout<< "to " << str.name() << endl;

    for
    (
        typename Triangulation::Finite_vertices_iterator vit =
            t.finite_vertices_begin();
        vit != t.finite_vertices_end();
        ++vit
    )
    {
        if (vit->type() >= startPointType && vit->type() <= endPointType)
        {
            meshTools::writeOBJ(str, topoint(vit->point()));
        }
    }
}


template<class Triangulation>
void Foam::DelaunayMeshTools::writeOBJ
(
    const fileName& fName,
    const Triangulation& t,
    const indexedVertexEnum::vertexType pointType
)
{
    writeOBJ(fName, t, pointType, pointType);
}


template<class Triangulation>
void Foam::DelaunayMeshTools::writeFixedPoints
(
    const fileName& fName,
    const Triangulation& t
)
{
    OFstream str(fName);

    Pout<< nl
        << "Writing fixed points to " << str.name() << endl;

    for
    (
        typename Triangulation::Finite_vertices_iterator vit =
            t.finite_vertices_begin();
        vit != t.finite_vertices_end();
        ++vit
    )
    {
        if (vit->fixed())
        {
            meshTools::writeOBJ(str, topoint(vit->point()));
        }
    }
}


template<class Triangulation>
void Foam::DelaunayMeshTools::writeBoundaryPoints
(
    const fileName& fName,
    const Triangulation& t
)
{
    OFstream str(fName);

    Pout<< nl
        << "Writing boundary points to " << str.name() << endl;

    for
    (
        typename Triangulation::Finite_vertices_iterator vit =
            t.finite_vertices_begin();
        vit != t.finite_vertices_end();
        ++vit
    )
    {
        if (!vit->internalPoint())
        {
            meshTools::writeOBJ(str, topoint(vit->point()));
        }
    }
}


template<class Triangulation>
void Foam::DelaunayMeshTools::writeProcessorInterface
(
    const fileName& fName,
    const Triangulation& t,
    const faceList& faces
)
{
    OFstream str(fName);

    pointField points(t.number_of_finite_cells(), point::max);

    for
    (
        typename Triangulation::Finite_cells_iterator cit =
            t.finite_cells_begin();
        cit != t.finite_cells_end();
        ++cit
    )
    {
        if (!cit->hasFarPoint() && !t.is_infinite(cit))
        {
            points[cit->cellIndex()] = cit->dual();
        }
    }

    meshTools::writeOBJ(str, faces, points);
}


template<class Triangulation>
void Foam::DelaunayMeshTools::writeInternalDelaunayVertices
(
    const fileName& instance,
    const Triangulation& t
)
{
    pointField internalDelaunayVertices(t.number_of_vertices());

    label vertI = 0;

    for
    (
        typename Triangulation::Finite_vertices_iterator vit =
            t.finite_vertices_begin();
        vit != t.finite_vertices_end();
        ++vit
    )
    {
        if (vit->internalPoint())
        {
            internalDelaunayVertices[vertI++] = topoint(vit->point());
        }
    }

    internalDelaunayVertices.setSize(vertI);

    pointIOField internalDVs
    (
        IOobject
        (
            "internalDelaunayVertices",
            instance,
            t.time(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        internalDelaunayVertices
    );

    Info<< nl
        << "Writing " << internalDVs.name()
        << " to " << internalDVs.instance()
        << endl;

    internalDVs.write();
}


template<class CellHandle>
void Foam::DelaunayMeshTools::drawDelaunayCell
(
    Ostream& os,
    const CellHandle& c,
    label offset
)
{
    // Supply offset as tet number
    offset *= 4;

    os  << "# cell index: " << label(c->cellIndex())
        << " INT_MIN = " << INT_MIN
        << endl;

    os  << "# circumradius "
        << mag(c->dual() - topoint(c->vertex(0)->point()))
        << endl;

    for (int i = 0; i < 4; i++)
    {
        os  << "# index / type / procIndex: "
            << label(c->vertex(i)->index()) << " "
            << label(c->vertex(i)->type()) << " "
            << label(c->vertex(i)->procIndex())
            <<
                (
                    CGAL::indexedVertexOps::uninitialised(c->vertex(i))
                  ? " # This vertex is uninitialised!"
                  : ""
                )
            << endl;

        meshTools::writeOBJ(os, topoint(c->vertex(i)->point()));
    }

    os  << "f " << 1 + offset << " " << 3 + offset << " " << 2 + offset << nl
        << "f " << 2 + offset << " " << 3 + offset << " " << 4 + offset << nl
        << "f " << 1 + offset << " " << 4 + offset << " " << 3 + offset << nl
        << "f " << 1 + offset << " " << 2 + offset << " " << 4 + offset << endl;

//    os  << "# cicumcentre " << endl;

//    meshTools::writeOBJ(os, c->dual());

//    os  << "l " << 1 + offset << " " << 5 + offset << endl;
}


template<class Triangulation>
Foam::tmp<Foam::pointField> Foam::DelaunayMeshTools::allPoints
(
    const Triangulation& t
)
{
    tmp<pointField> tpts(new pointField(t.vertexCount(), point::max));
    pointField& pts = tpts.ref();

    for
    (
        typename Triangulation::Finite_vertices_iterator vit =
            t.finite_vertices_begin();
        vit != t.finite_vertices_end();
        ++vit
    )
    {
        if (vit->internalOrBoundaryPoint() && !vit->referred())
        {
            pts[vit->index()] = topoint(vit->point());
        }
    }

    return tpts;
}


// ************************************************************************* //
