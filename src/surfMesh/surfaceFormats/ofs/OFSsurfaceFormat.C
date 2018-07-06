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

#include "OFSsurfaceFormat.H"
#include "IFstream.H"
#include "IStringStream.H"
#include "ListOps.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Face>
Foam::fileFormats::OFSsurfaceFormat<Face>::OFSsurfaceFormat
(
    const fileName& filename
)
{
    read(filename);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Face>
bool Foam::fileFormats::OFSsurfaceFormat<Face>::read
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

    // read surfZones:
    is >> this->storedZones();

    // read points:
    is >> this->storedPoints();

    // must triangulate?
    if (MeshedSurface<Face>::isTri())
    {
        // read faces as 'face' and transcribe to 'triFace'
        List<face> faceLst(is);

        MeshedSurface<face> surf
        (
            xferMove(this->storedPoints()),
            xferMove(faceLst),
            xferMove(this->storedZones())
        );

        this->transcribe(surf);
    }
    else
    {
        // read faces directly
        is >> this->storedFaces();
    }

    return true;
}


template<class Face>
bool Foam::fileFormats::OFSsurfaceFormat<Face>::read
(
    Istream& is,
    pointField& pointLst,
    List<Face>& faceLst,
    List<surfZone>& zoneLst
)
{
    if (!is.good())
    {
        FatalErrorInFunction
            << "read error "
            << exit(FatalError);
    }

    // read surfZones:
    is >> zoneLst;

    // read points:
    is >> pointLst;

    // must triangulate?
    if (MeshedSurface<Face>::isTri())
    {
        // read faces as 'face' and transcribe to 'triFace'
        List<face> origFaces(is);

        MeshedSurface<face> origSurf
        (
            xferMove(pointLst),
            xferMove(origFaces),
            xferMove(zoneLst)
        );

        MeshedSurface<Face> surf;
        surf.transcribe(origSurf);
    }
    else
    {
        // read faces directly
        is >> faceLst;
    }

    return true;
}


template<class Face>
bool Foam::fileFormats::OFSsurfaceFormat<Face>::read
(
    Istream& is,
    MeshedSurface<Face>& surf
)
{
    surf.clear();

    if (!is.good())
    {
        FatalErrorInFunction
            << "read error "
            << exit(FatalError);
    }

    pointField pointLst;
    List<Face> faceLst;
    List<surfZone> zoneLst;

    read(is, pointLst, faceLst, zoneLst);

    surf.reset
    (
        xferMove(pointLst),
        xferMove(faceLst),
        xferMove(zoneLst)
    );

    return true;
}


template<class Face>
bool Foam::fileFormats::OFSsurfaceFormat<Face>::read
(
    Istream& is,
    UnsortedMeshedSurface<Face>& surf
)
{
    surf.clear();
    MeshedSurface<Face> origSurf(is);
    surf.transfer(origSurf);

    return true;
}



template<class Face>
void Foam::fileFormats::OFSsurfaceFormat<Face>::write
(
    const fileName& filename,
    const MeshedSurfaceProxy<Face>& surf
)
{
    const List<Face>&  faceLst = surf.faces();
    const List<label>& faceMap = surf.faceMap();

    OFstream os(filename);
    if (!os.good())
    {
        FatalErrorInFunction
            << "Cannot open file for writing " << filename
            << exit(FatalError);
    }


    OFSsurfaceFormatCore::writeHeader(os, surf.points(), surf.surfZones());

    const List<surfZone>& zones = surf.surfZones();
    const bool useFaceMap = (surf.useFaceMap() && zones.size() > 1);

    if (useFaceMap)
    {
        os  << "\n// faces:"  << nl
            << faceLst.size() << token::BEGIN_LIST << nl;

        label facei = 0;
        forAll(zones, zoneI)
        {
            // Print all faces belonging to this zone
            const surfZone& zone = zones[zoneI];

            forAll(zone, localFacei)
            {
                os << faceLst[faceMap[facei++]] << nl;
            }
        }
        os << token::END_LIST << nl;
    }
    else
    {
        os  << "\n// faces:"  << nl << faceLst << nl;
    }

    IOobject::writeDivider(os);

    // Check state of Ostream
    os.check("OFSsurfaceFormat<Face>::write(Ostream&)");
}


// ************************************************************************* //
