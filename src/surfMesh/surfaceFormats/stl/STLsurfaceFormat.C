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

#include "STLsurfaceFormat.H"
#include "ListOps.H"
#include "triPointRef.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Face>
inline void Foam::fileFormats::STLsurfaceFormat<Face>::writeShell
(
    Ostream& os,
    const pointField& pointLst,
    const Face& f
)
{
    // calculate the normal ourselves, for flexibility and speed
    const vector norm = triPointRef
    (
        pointLst[f[0]],
        pointLst[f[1]],
        pointLst[f[2]]
    ).normal();

    // simple triangulation about f[0].
    // better triangulation should have been done before
    const point& p0 = pointLst[f[0]];
    for (label fp1 = 1; fp1 < f.size() - 1; ++fp1)
    {
        label fp2 = f.fcIndex(fp1);

        const point& p1 = pointLst[f[fp1]];
        const point& p2 = pointLst[f[fp2]];

        // write STL triangle
        os  << " facet normal "
            << norm.x() << ' ' << norm.y() << ' ' << norm.z() << nl
            << "  outer loop\n"
            << "   vertex " << p0.x() << ' ' << p0.y() << ' ' << p0.z() << nl
            << "   vertex " << p1.x() << ' ' << p1.y() << ' ' << p1.z() << nl
            << "   vertex " << p2.x() << ' ' << p2.y() << ' ' << p2.z() << nl
            << "  endloop\n"
            << " endfacet" << endl;
    }
}


template<class Face>
inline void Foam::fileFormats::STLsurfaceFormat<Face>::writeShell
(
    ostream& os,
    const pointField& pointLst,
    const Face& f,
    const label zoneI
)
{
    // calculate the normal ourselves, for flexibility and speed
    const vector norm = triPointRef
    (
        pointLst[f[0]],
        pointLst[f[1]],
        pointLst[f[2]]
    ).normal();

    // simple triangulation about f[0].
    // better triangulation should have been done before
    const point& p0 = pointLst[f[0]];
    for (label fp1 = 1; fp1 < f.size() - 1; ++fp1)
    {
        label fp2 = f.fcIndex(fp1);

        STLtriangle stlTri
        (
            norm,
            p0,
            pointLst[f[fp1]],
            pointLst[f[fp2]],
            zoneI
        );

        stlTri.write(os);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Face>
Foam::fileFormats::STLsurfaceFormat<Face>::STLsurfaceFormat
(
    const fileName& filename
)
{
    read(filename);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Face>
bool Foam::fileFormats::STLsurfaceFormat<Face>::read
(
    const fileName& filename
)
{
    this->clear();

    // read in the values
    STLsurfaceFormatCore reader(filename);

    // transfer points
    this->storedPoints().transfer(reader.points());

    // retrieve the original zone information
    List<word>  names(reader.names().xfer());
    List<label> sizes(reader.sizes().xfer());
    List<label> zoneIds(reader.zoneIds().xfer());

    // generate the (sorted) faces
    List<Face> faceLst(zoneIds.size());

    if (reader.sorted())
    {
        // already sorted - generate directly
        forAll(faceLst, facei)
        {
            const label startPt = 3*facei;
            faceLst[facei] = triFace(startPt, startPt+1, startPt+2);
        }
    }
    else
    {
        // unsorted - determine the sorted order:
        // avoid SortableList since we discard the main list anyhow
        List<label> faceMap;
        sortedOrder(zoneIds, faceMap);

        // generate sorted faces
        forAll(faceMap, facei)
        {
            const label startPt = 3*faceMap[facei];
            faceLst[facei] = triFace(startPt, startPt+1, startPt+2);
        }
    }
    zoneIds.clear();

    // transfer:
    this->storedFaces().transfer(faceLst);

    if (names.size())
    {
        this->addZones(sizes, names);
    }
    else
    {
        this->addZones(sizes);
    }

    this->stitchFaces(small);
    return true;
}



template<class Face>
void Foam::fileFormats::STLsurfaceFormat<Face>::writeAscii
(
    const fileName& filename,
    const MeshedSurfaceProxy<Face>& surf
)
{
    OFstream os(filename);
    if (!os.good())
    {
        FatalErrorInFunction
            << "Cannot open file for writing " << filename
            << exit(FatalError);
    }

    const pointField& pointLst = surf.points();
    const List<Face>&  faceLst = surf.faces();
    const List<label>& faceMap = surf.faceMap();

    const List<surfZone>& zones =
    (
        surf.surfZones().empty()
      ? surfaceFormatsCore::oneZone(faceLst)
      : surf.surfZones()
    );

    const bool useFaceMap = (surf.useFaceMap() && zones.size() > 1);

    label faceIndex = 0;
    forAll(zones, zoneI)
    {
        // Print all faces belonging to this zone
        const surfZone& zone = zones[zoneI];

        os << "solid " << zone.name() << nl;

        if (useFaceMap)
        {
            forAll(zone, localFacei)
            {
                const label facei = faceMap[faceIndex++];
                writeShell(os, pointLst, faceLst[facei]);
            }
        }
        else
        {
            forAll(zone, localFacei)
            {
                writeShell(os, pointLst, faceLst[faceIndex++]);
            }
        }
        os << "endsolid " << zone.name() << endl;
    }
}


template<class Face>
void Foam::fileFormats::STLsurfaceFormat<Face>::writeBinary
(
    const fileName& filename,
    const MeshedSurfaceProxy<Face>& surf
)
{
    std::ofstream os(filename.c_str(), std::ios::binary);
    if (!os.good())
    {
        FatalErrorInFunction
            << "Cannot open file for writing " << filename
            << exit(FatalError);
    }


    const pointField& pointLst = surf.points();
    const List<Face>&  faceLst = surf.faces();
    const List<label>& faceMap = surf.faceMap();

    const List<surfZone>& zones =
    (
        surf.surfZones().size() > 1
      ? surf.surfZones()
      : surfaceFormatsCore::oneZone(faceLst)
    );

    const bool useFaceMap = (surf.useFaceMap() && zones.size() > 1);


    unsigned int nTris = 0;
    if (MeshedSurface<Face>::isTri())
    {
        nTris = faceLst.size();
    }
    else
    {
        // count triangles for on-the-fly triangulation
        forAll(faceLst, facei)
        {
            nTris += faceLst[facei].size() - 2;
        }
    }

    // Write the STL header
    STLsurfaceFormatCore::writeHeaderBINARY(os, nTris);

    label faceIndex = 0;
    forAll(zones, zoneI)
    {
        const surfZone& zone = zones[zoneI];

        if (useFaceMap)
        {
            forAll(zone, localFacei)
            {
                writeShell
                (
                    os,
                    pointLst,
                    faceLst[faceMap[faceIndex++]],
                    zoneI
                );
            }
        }
        else
        {
            forAll(zone, localFacei)
            {
                writeShell
                (
                    os,
                    pointLst,
                    faceLst[faceIndex++],
                    zoneI
                );
            }
        }
    }
}


template<class Face>
void Foam::fileFormats::STLsurfaceFormat<Face>::writeAscii
(
    const fileName& filename,
    const UnsortedMeshedSurface<Face>& surf
)
{
    OFstream os(filename);
    if (!os.good())
    {
        FatalErrorInFunction
            << "Cannot open file for writing " << filename
            << exit(FatalError);
    }

    // a single zone - we can skip sorting
    if (surf.zoneToc().size() == 1)
    {
        const pointField& pointLst = surf.points();
        const List<Face>& faceLst  = surf.faces();

        os << "solid " << surf.zoneToc()[0].name() << endl;
        forAll(faceLst, facei)
        {
            writeShell(os, pointLst, faceLst[facei]);
        }
        os << "endsolid " << surf.zoneToc()[0].name() << endl;
    }
   else
   {
       labelList faceMap;
       List<surfZone> zoneLst = surf.sortedZones(faceMap);

       writeAscii
       (
           filename,
           MeshedSurfaceProxy<Face>
           (
               surf.points(),
               surf.faces(),
               zoneLst,
               faceMap
           )
       );
   }
}


template<class Face>
void Foam::fileFormats::STLsurfaceFormat<Face>::writeBinary
(
    const fileName& filename,
    const UnsortedMeshedSurface<Face>& surf
)
{
    std::ofstream os(filename.c_str(), std::ios::binary);
    if (!os.good())
    {
        FatalErrorInFunction
            << "Cannot open file for writing " << filename
            << exit(FatalError);
    }

    const pointField&  pointLst = surf.points();
    const List<Face>&  faceLst  = surf.faces();
    const List<label>& zoneIds  = surf.zoneIds();

    unsigned int nTris = 0;
    if (MeshedSurface<Face>::isTri())
    {
        nTris = faceLst.size();
    }
    else
    {
        // count triangles for on-the-fly triangulation
        forAll(faceLst, facei)
        {
            nTris += faceLst[facei].size() - 2;
        }
    }

    // Write the STL header
    STLsurfaceFormatCore::writeHeaderBINARY(os, nTris);

    // always write unsorted
    forAll(faceLst, facei)
    {
        writeShell
        (
            os,
            pointLst,
            faceLst[facei],
            zoneIds[facei]
        );
    }
}


template<class Face>
void Foam::fileFormats::STLsurfaceFormat<Face>::write
(
    const fileName& filename,
    const MeshedSurfaceProxy<Face>& surf
)
{
    const word ext = filename.ext();

    // handle 'stlb' as binary directly
    if (ext == "stlb")
    {
        writeBinary(filename, surf);
    }
    else
    {
        writeAscii(filename, surf);
    }
}


template<class Face>
void Foam::fileFormats::STLsurfaceFormat<Face>::write
(
    const fileName& filename,
    const UnsortedMeshedSurface<Face>& surf
)
{
    word ext = filename.ext();

    // handle 'stlb' as binary directly
    if (ext == "stlb")
    {
        writeBinary(filename, surf);
    }
    else
    {
        writeAscii(filename, surf);
    }
}


// ************************************************************************* //
