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

#include "VTKsurfaceFormat.H"
#include "vtkUnstructuredReader.H"
#include "scalarIOField.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Face>
void Foam::fileFormats::VTKsurfaceFormat<Face>::writeHeaderPolygons
(
    Ostream& os,
    const UList<Face>& faceLst
)
{
    label nNodes = 0;

    forAll(faceLst, facei)
    {
        nNodes += faceLst[facei].size();
    }

    os  << nl
        << "POLYGONS " << faceLst.size() << ' '
        << faceLst.size() + nNodes << nl;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Face>
Foam::fileFormats::VTKsurfaceFormat<Face>::VTKsurfaceFormat
(
    const fileName& filename
)
{
    read(filename);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Face>
bool Foam::fileFormats::VTKsurfaceFormat<Face>::read
(
    const fileName& filename
)
{
    const bool mustTriangulate = this->isTri();
    this->clear();

    IFstream is(filename);
    if (!is.good())
    {
        FatalErrorInFunction
            << "Cannot read file " << filename
            << exit(FatalError);
    }

    // assume that the groups are not intermixed
    bool sorted = true;


    // Construct dummy time so we have something to create an objectRegistry
    // from
    Time dummyTime
    (
        "dummyRoot",
        "dummyCase",
        "system",
        "constant",
        false           // enableFunctionObjects
    );

    // Make dummy object registry
    objectRegistry obr
    (
        IOobject
        (
            "dummy",
            dummyTime,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        )
    );

    // Read all
    vtkUnstructuredReader reader(obr, is);
    const faceList& faces = reader.faces();

    // Assume all faces in zone0 unless a region field is present
    labelList zones(faces.size(), 0);
    if (reader.cellData().foundObject<scalarIOField>("region"))
    {
        const scalarIOField& region =
            reader.cellData().lookupObject<scalarIOField>
            (
                "region"
            );
        forAll(region, i)
        {
            zones[i] = label(region[i]);
        }
    }
    else if (reader.cellData().foundObject<scalarIOField>("STLSolidLabeling"))
    {
        const scalarIOField& region =
            reader.cellData().lookupObject<scalarIOField>
            (
                "STLSolidLabeling"
            );
        forAll(region, i)
        {
            zones[i] = label(region[i]);
        }
    }

    // Create zone names
    const label nZones = max(zones)+1;
    wordList zoneNames(nZones);
    forAll(zoneNames, i)
    {
        zoneNames[i] = "zone" + Foam::name(i);
    }


    // See if needs triangulation
    label nTri = 0;
    if (mustTriangulate)
    {
        forAll(faces, facei)
        {
            nTri += faces[facei].size()-2;
        }
    }

    if (nTri > 0)
    {
        DynamicList<Face> dynFaces(nTri);
        DynamicList<label> dynZones(nTri);
        forAll(faces, facei)
        {
            const face& f = faces[facei];
            for (label fp1 = 1; fp1 < f.size() - 1; fp1++)
            {
                label fp2 = f.fcIndex(fp1);

                dynFaces.append(triFace(f[0], f[fp1], f[fp2]));
                dynZones.append(zones[facei]);
            }
        }

        // Count
        labelList zoneSizes(nZones, 0);
        forAll(dynZones, triI)
        {
            zoneSizes[dynZones[triI]]++;
        }

        this->sortFacesAndStore(dynFaces.xfer(), dynZones.xfer(), sorted);

        // add zones, culling empty ones
        this->addZones(zoneSizes, zoneNames, true);
    }
    else
    {
        DynamicList<Face> dynFaces(faces.size());
        forAll(faces, facei)
        {
            const face& f = faces[facei];
            dynFaces.append(Face(f));
        }

        // Count
        labelList zoneSizes(nZones, 0);
        forAll(zones, facei)
        {
            zoneSizes[zones[facei]]++;
        }

        this->sortFacesAndStore(dynFaces.xfer(), zones.xfer(), sorted);

        // add zones, culling empty ones
        this->addZones(zoneSizes, zoneNames, true);
    }

    // transfer to normal lists
    this->storedPoints().transfer(reader.points());

    return true;
}


template<class Face>
void Foam::fileFormats::VTKsurfaceFormat<Face>::write
(
    const fileName& filename,
    const MeshedSurfaceProxy<Face>& surf
)
{
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

    OFstream os(filename);
    if (!os.good())
    {
        FatalErrorInFunction
            << "Cannot open file for writing " << filename
            << exit(FatalError);
    }


    writeHeader(os, pointLst);
    writeHeaderPolygons(os, faceLst);

    label faceIndex = 0;
    forAll(zones, zoneI)
    {
        const surfZone& zone = zones[zoneI];

        if (useFaceMap)
        {
            forAll(zone, localFacei)
            {
                const Face& f = faceLst[faceMap[faceIndex++]];

                os << f.size();
                forAll(f, fp)
                {
                    os << ' ' << f[fp];
                }
                os << ' ' << nl;
            }
        }
        else
        {
            forAll(zone, localFacei)
            {
                const Face& f = faceLst[faceIndex++];

                os << f.size();
                forAll(f, fp)
                {
                    os << ' ' << f[fp];
                }
                os << ' ' << nl;
            }
        }
    }

    writeTail(os, zones);
}


template<class Face>
void Foam::fileFormats::VTKsurfaceFormat<Face>::write
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


    const List<Face>& faceLst = surf.faces();

    writeHeader(os, surf.points());
    writeHeaderPolygons(os, faceLst);

    forAll(faceLst, facei)
    {
        const Face& f = faceLst[facei];

        os << f.size();
        forAll(f, fp)
        {
            os << ' ' << f[fp];
        }
        os << ' ' << nl;
    }

    writeTail(os, surf.zoneIds());
}


// ************************************************************************* //
