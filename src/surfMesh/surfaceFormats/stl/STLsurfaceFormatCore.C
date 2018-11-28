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

#include "STLsurfaceFormatCore.H"
#include "gzstream.h"
#include "OSspecific.H"
#include "Map.H"
#include "IFstream.H"
#include "Ostream.H"

#undef DEBUG_STLBINARY

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// check binary by getting the header and number of facets
// this seems to work better than the old token-based method
// - some programs (eg, pro-STAR) have 'solid' as the first word in
//   the binary header.
// - using wordToken can cause an abort if non-word (binary) content
//   is detected ... this is not exactly what we want.
int Foam::fileFormats::STLsurfaceFormatCore::detectBINARY
(
    const fileName& filename
)
{
    off_t dataFileSize = Foam::fileSize(filename);

    IFstream str(filename, IOstream::BINARY);
    istream& is = str().stdStream();

    // Read the STL header
    char header[headerSize];
    is.read(header, headerSize);

    // Check that stream is OK, if not this may be an ASCII file
    if (!is.good())
    {
        return 0;
    }

    // Read the number of triangles in the STl file
    // (note: read as int so we can check whether >2^31)
    int nTris;
    is.read(reinterpret_cast<char*>(&nTris), sizeof(unsigned int));

    // Check that stream is OK and number of triangles is positive,
    // if not this may be an ASCII file
    //
    // Also compare the file size with that expected from the number of tris
    // If the comparison is not sensible then it may be an ASCII file
    if
    (
        !is
     || nTris < 0
     || nTris < (dataFileSize - headerSize)/50
     || nTris > (dataFileSize - headerSize)/25
    )
    {
        return 0;
    }

    // looks like it might be BINARY, return number of triangles
    return nTris;
}


bool Foam::fileFormats::STLsurfaceFormatCore::readBINARY
(
    istream& is,
    const off_t dataFileSize
)
{
    sorted_ = true;

    // Read the STL header
    char header[headerSize];
    is.read(header, headerSize);

    // Check that stream is OK, if not this may be an ASCII file
    if (!is.good())
    {
        FatalErrorInFunction
            << "problem reading header, perhaps file is not binary "
            << exit(FatalError);
    }

    // Read the number of triangles in the STl file
    // (note: read as int so we can check whether >2^31)
    int nTris;
    is.read(reinterpret_cast<char*>(&nTris), sizeof(unsigned int));

    // Check that stream is OK and number of triangles is positive,
    // if not this maybe an ASCII file
    //
    // Also compare the file size with that expected from the number of tris
    // If the comparison is not sensible then it may be an ASCII file
    if
    (
        !is
     || nTris < 0
     || nTris < int(dataFileSize - headerSize)/50
     || nTris > int(dataFileSize - headerSize)/25
    )
    {
        FatalErrorInFunction
            << "problem reading number of triangles, perhaps file is not binary"
            << exit(FatalError);
    }

#ifdef DEBUG_STLBINARY
    Info<< "# " << nTris << " facets" << endl;
    label prevZone = -1;
#endif

    points_.setSize(3*nTris);
    zoneIds_.setSize(nTris);

    Map<label> lookup;
    DynamicList<label> dynSizes;

    label ptI = 0;
    label zoneI = -1;
    forAll(zoneIds_, facei)
    {
        // Read an STL triangle
        STLtriangle stlTri(is);

        // transcribe the vertices of the STL triangle -> points
        points_[ptI++] = stlTri.a();
        points_[ptI++] = stlTri.b();
        points_[ptI++] = stlTri.c();

        // interpret stl attribute as a zone
        const label origId = stlTri.attrib();

        Map<label>::const_iterator fnd = lookup.find(origId);
        if (fnd != lookup.end())
        {
            if (zoneI != fnd())
            {
                // group appeared out of order
                sorted_ = false;
            }
            zoneI = fnd();
        }
        else
        {
            zoneI = dynSizes.size();
            lookup.insert(origId, zoneI);
            dynSizes.append(0);
        }

        zoneIds_[facei] = zoneI;
        dynSizes[zoneI]++;

#ifdef DEBUG_STLBINARY
        if (prevZone != zoneI)
        {
            if (prevZone != -1)
            {
                Info<< "endsolid zone" << prevZone << nl;
            }
            prevZone = zoneI;

            Info<< "solid zone" << prevZone << nl;
        }

        Info<< " facet normal " << stlTri.normal() << nl
            << "  outer loop" << nl
            << "   vertex " << stlTri.a() << nl
            << "   vertex " << stlTri.b() << nl
            << "   vertex " << stlTri.c() << nl
            << "  outer loop" << nl
            << " endfacet" << endl;
#endif
    }

    names_.clear();
    sizes_.transfer(dynSizes);

    return true;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fileFormats::STLsurfaceFormatCore::STLsurfaceFormatCore
(
    const fileName& filename
)
:
    sorted_(true),
    points_(0),
    zoneIds_(0),
    names_(0),
    sizes_(0)
{
    off_t dataFileSize = Foam::fileSize(filename);

    // auto-detect ascii/binary
    if (detectBINARY(filename))
    {
        readBINARY
        (
            IFstream(filename, IOstream::BINARY)().stdStream(),
            dataFileSize
        );
    }
    else
    {
        readASCII
        (
            IFstream(filename)().stdStream(),
            dataFileSize
        );
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fileFormats::STLsurfaceFormatCore::~STLsurfaceFormatCore()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fileFormats::STLsurfaceFormatCore::writeHeaderBINARY
(
    ostream& os,
    unsigned int nTris
)
{
    // STL header with extra information about nTris
    char header[headerSize];
    sprintf(header, "STL binary file %u facets", nTris);

    // avoid trailing junk
    for (size_t i = strlen(header); i < headerSize; ++i)
    {
        header[i] = 0;
    }

    os.write(header, headerSize);
    os.write(reinterpret_cast<char*>(&nTris), sizeof(unsigned int));

}


// ************************************************************************* //
