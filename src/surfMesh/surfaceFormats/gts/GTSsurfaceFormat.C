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

#include "GTSsurfaceFormat.H"
#include "surfaceFormatsCore.H"
#include "clock.H"
#include "IFstream.H"
#include "IStringStream.H"
#include "Ostream.H"
#include "OFstream.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Face>
Foam::fileFormats::GTSsurfaceFormat<Face>::GTSsurfaceFormat
(
    const fileName& filename
)
{
    read(filename);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Face>
bool Foam::fileFormats::GTSsurfaceFormat<Face>::read
(
    const fileName& filename
)
{
    this->clear();

    IFstream is(filename);
    if (!is.good())
    {
        FatalErrorInFunction
            << "Cannot read file " << filename
            << exit(FatalError);
    }

    // Read header
    string line = this->getLineNoComment(is);

    label nPoints, nEdges, nElems;
    {
        IStringStream lineStream(line);
        lineStream
            >> nPoints
            >> nEdges
            >> nElems;
    }


    // write directly into the lists:
    pointField&  pointLst = this->storedPoints();
    List<Face>&  faceLst  = this->storedFaces();
    List<label>& zoneIds  = this->storedZoneIds();

    pointLst.setSize(nPoints);
    faceLst.setSize(nElems);
    zoneIds.setSize(nElems);

    // Read points
    forAll(pointLst, pointi)
    {
        scalar x, y, z;
        line = this->getLineNoComment(is);
        {
            IStringStream lineStream(line);
            lineStream
                >> x >> y >> z;
        }

        pointLst[pointi] = point(x, y, z);
    }

    // Read edges (Foam indexing)
    edgeList edges(nEdges);
    forAll(edges, edgei)
    {
        label beg, end;
        line = this->getLineNoComment(is);
        {
            IStringStream lineStream(line);
            lineStream
                >> beg >> end;
        }
        edges[edgei] = edge(beg - 1, end - 1);
    }


    // Read triangles. Convert references to edges into pointlabels
    label maxZone = 0;
    forAll(faceLst, facei)
    {
        label e0Label, e1Label, e2Label;
        label zoneI = 0;

        line = this->getLineNoComment(is);
        {
            IStringStream lineStream(line);
            lineStream
                >> e0Label >> e1Label >> e2Label;

            // Optional zone number: read first, then check state on stream
            if (lineStream)
            {
                label num;
                lineStream >> num;
                if (!lineStream.bad())
                {
                    zoneI = num;
                    if (maxZone < zoneI)
                    {
                        maxZone = zoneI;
                    }
                }
            }
        }

        // Determine ordering of edges e0, e1
        //  common: common vertex, shared by e0 and e1
        //  e0Far:  vertex on e0 which is not common
        //  e1Far:  vertex on e1 which is not common
        const edge& e0 = edges[e0Label - 1];
        const edge& e1 = edges[e1Label - 1];
        const edge& e2 = edges[e2Label - 1];

        label common01 = e0.commonVertex(e1);
        if (common01 == -1)
        {
            FatalErrorInFunction
                << "Edges 0 and 1 of triangle " << facei
                << " do not share a point.\n"
                << "    edge0:" << e0 << nl
                << "    edge1:" << e1
                << exit(FatalError);
        }

        label e0Far = e0.otherVertex(common01);
        label e1Far = e1.otherVertex(common01);

        label common12 = e1.commonVertex(e2);
        if (common12 == -1)
        {
            FatalErrorInFunction
                << "Edges 1 and 2 of triangle " << facei
                << " do not share a point.\n"
                << "    edge1:" << e1 << nl
                << "    edge2:" << e2
                << exit(FatalError);
        }
        label e2Far = e2.otherVertex(common12);

        // Does edge2 sit between edge1 and 0?
        if (common12 != e1Far || e2Far != e0Far)
        {
            FatalErrorInFunction
                << "Edges of triangle " << facei
                << " reference more than three points.\n"
                << "    edge0:" << e0 << nl
                << "    edge1:" << e1 << nl
                << "    edge2:" << e2 << nl
                << exit(FatalError);
        }

        faceLst[facei] = triFace(e0Far, common01, e1Far);
        zoneIds[facei] = zoneI;
    }


    List<surfZoneIdentifier> newZones(maxZone+1);
    forAll(newZones, zoneI)
    {
        newZones[zoneI] = surfZoneIdentifier
        (
            "zone" + ::Foam::name(zoneI),
            zoneI
        );
    }

    this->storedZoneToc().transfer(newZones);

    return true;
}


template<class Face>
void Foam::fileFormats::GTSsurfaceFormat<Face>::write
(
    const fileName& filename,
    const MeshedSurface<Face>& surf
)
{
    const pointField& pointLst = surf.points();
    const List<Face>& faceLst  = surf.faces();

    const List<surfZone>& zones =
    (
        surf.surfZones().size()
      ? surf.surfZones()
      : surfaceFormatsCore::oneZone(faceLst)
    );

    // check if output triangulation would be required
    // It is too annoying to triangulate on-the-fly
    // just issue a warning and get out
    if (!MeshedSurface<Face>::isTri())
    {
        label nNonTris = 0;
        forAll(faceLst, facei)
        {
            if (faceLst[facei].size() != 3)
            {
                ++nNonTris;
            }
        }

        if (nNonTris)
        {
            FatalErrorInFunction
                << "Surface has " << nNonTris << "/" << faceLst.size()
                << " non-triangulated faces - not writing!" << endl;
            return;
        }
    }


    OFstream os(filename);
    if (!os.good())
    {
        FatalErrorInFunction
            << "Cannot open file for writing " << filename
            << exit(FatalError);
    }


    // Write header, print zone names as comment
    os  << "# GTS file" << nl
        << "# Zones:" << nl;

    forAll(zones, zoneI)
    {
        os  << "#     " << zoneI << "    "
            << zones[zoneI].name() << nl;
    }
    os  << "#" << nl;

    os  << "# nPoints  nEdges  nTriangles" << nl
        << pointLst.size() << ' ' << surf.nEdges() << ' '
        << surf.size() << endl;


    // Write vertex coords
    forAll(pointLst, pointi)
    {
        const point& pt = pointLst[pointi];

        os  << pt.x() << ' ' << pt.y() << ' ' << pt.z() << nl;
    }


    // Write edges.
    // Note: edges are in local point labels so convert
    const edgeList& es = surf.edges();
    const labelList& meshPts = surf.meshPoints();

    forAll(es, edgei)
    {
        os  << meshPts[es[edgei].start()] + 1 << ' '
            << meshPts[es[edgei].end()] + 1 << endl;
    }

    // Write faces in terms of edges.
    const labelListList& faceEs = surf.faceEdges();

    label faceIndex = 0;
    forAll(zones, zoneI)
    {
        const surfZone& zone = zones[zoneI];

        forAll(zone, localFacei)
        {
            const labelList& fEdges = faceEs[faceIndex++];

            os  << fEdges[0] + 1 << ' '
                << fEdges[1] + 1 << ' '
                << fEdges[2] + 1 << ' '
                << zoneI << endl;
        }
    }
}


template<class Face>
void Foam::fileFormats::GTSsurfaceFormat<Face>::write
(
    const fileName& filename,
    const UnsortedMeshedSurface<Face>& surf
)
{
    const pointField& pointLst   = surf.points();
    const List<Face>& faceLst    = surf.faces();
    const List<label>& zoneIds   = surf.zoneIds();
    const List<surfZoneIdentifier>& zoneToc = surf.zoneToc();

    // check if output triangulation would be required
    // It is too annoying to triangulate on-the-fly
    // just issue a warning and get out
    if (!MeshedSurface<Face>::isTri())
    {
        label nNonTris = 0;
        forAll(faceLst, facei)
        {
            if (faceLst[facei].size() != 3)
            {
                ++nNonTris;
            }
        }

        if (nNonTris)
        {
            FatalErrorInFunction
                << "Surface has " << nNonTris << "/" << faceLst.size()
                << " non-triangulated faces - not writing!" << endl;
            return;
        }
    }


    OFstream os(filename);
    if (!os.good())
    {
        FatalErrorInFunction
            << "Cannot open file for writing " << filename
            << exit(FatalError);
    }


    // Write header, print zone names as comment
    os  << "# GTS file" << nl
        << "# Zones:" << nl;

    forAll(zoneToc, zoneI)
    {
        os  << "#     " << zoneI << "    "
            << zoneToc[zoneI].name() << nl;
    }
    os  << "#" << endl;


    os  << "# nPoints  nEdges  nTriangles" << nl
        << pointLst.size() << ' ' << surf.nEdges() << ' '
        << surf.size() << endl;


    // Write vertex coords
    forAll(pointLst, pointi)
    {
        os  << pointLst[pointi].x() << ' '
            << pointLst[pointi].y() << ' '
            << pointLst[pointi].z() << endl;
    }


    // Write edges.
    // Note: edges are in local point labels so convert
    const edgeList& es = surf.edges();
    const labelList& meshPts = surf.meshPoints();

    forAll(es, edgeI)
    {
        os  << meshPts[es[edgeI].start()] + 1 << ' '
            << meshPts[es[edgeI].end()] + 1 << endl;
    }


    // Write faces in terms of edges.
    const labelListList& faceEs = surf.faceEdges();

    forAll(faceLst, facei)
    {
        const labelList& fEdges = faceEs[facei];

        os  << fEdges[0] + 1 << ' '
            << fEdges[1] + 1 << ' '
            << fEdges[2] + 1 << ' '
            << zoneIds[facei] << endl;
    }
}


// ************************************************************************* //
