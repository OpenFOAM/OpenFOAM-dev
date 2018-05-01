/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
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

#include "edgeMeshFormat.H"
#include "IOobject.H"
#include "IFstream.H"
#include "clock.H"
#include "Time.H"
#include "featureEdgeMesh.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fileFormats::edgeMeshFormat::edgeMeshFormat
(
    const fileName& filename
)
{
    read(filename);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::fileFormats::edgeMeshFormat::read
(
    const fileName& filename
)
{
    clear();

    fileName dir = filename.path();
    fileName caseName = dir.name();
    fileName rootPath = dir.path();

    // Construct dummy time to use as an objectRegistry
    Time dummyTime
    (
        ".",        // rootPath,
        ".",        // caseName,
        "system",   // systemName,
        "constant", // constantName,
        false       // enableFunctionObjects
    );

    // Construct IOobject to re-use the headerOk & readHeader
    // (so we can read ascii and binary)
    IOobject io
    (
        filename,
        dummyTime,
        IOobject::NO_READ,
        IOobject::NO_WRITE,
        false
    );

    if (!io.typeHeaderOk<featureEdgeMesh>(false))
    {
        FatalErrorInFunction
            << "Cannot read file " << filename
            << exit(FatalError);
    }

    const fileName fName(typeFilePath<featureEdgeMesh>(io));

    autoPtr<IFstream> isPtr(new IFstream(fName));
    bool ok = false;
    if (isPtr().good())
    {
        Istream& is = isPtr();
        ok = io.readHeader(is);

        if (ok)
        {
            ok = read(is, this->storedPoints(), this->storedEdges());
        }
    }

    return ok;
}


bool Foam::fileFormats::edgeMeshFormat::read
(
    Istream& is,
    pointField& pointLst,
    edgeList& edgeLst
)
{
    if (!is.good())
    {
        FatalErrorInFunction
            << "read error "
            << exit(FatalError);
    }

    // read points:
    is  >> pointLst;

    // read edges:
    is  >> edgeLst;

    return true;
}


Foam::Ostream& Foam::fileFormats::edgeMeshFormat::write
(
    Ostream& os,
    const pointField& pointLst,
    const edgeList& edgeLst
)
{
    if (!os.good())
    {
        FatalErrorInFunction
            << "bad output stream " << os.name()
            << exit(FatalError);
    }

    os  << "\n// points:" << nl << pointLst << nl
        << "\n// edges:" << nl << edgeLst << nl;

    IOobject::writeDivider(os);

    // Check state of Ostream
    os.check
    (
        "edgeMeshFormat::write"
        "(Ostream&, const pointField&, const edgeList&)"
    );

    return os;
}


void Foam::fileFormats::edgeMeshFormat::write
(
    const fileName& filename,
    const edgeMesh& mesh
)
{
    // Construct dummy time to use as an objectRegistry
    Time dummyTime
    (
        ".",        // rootPath,
        ".",        // caseName,
        "system",   // systemName,
        "constant", // constantName,
        false       // enableFunctionObjects
    );

    // Construct IOobject to re-use the writeHeader
    IOobject io
    (
        filename,
        dummyTime,
        IOobject::NO_READ,
        IOobject::NO_WRITE,
        false
    );
    io.note() = "written " + clock::dateTime();

    // Note: always write ascii
    autoPtr<OFstream> osPtr(new OFstream(filename));

    if (!osPtr().good())
    {
        FatalIOErrorInFunction
        (
            osPtr()
        )   << "Cannot open file for writing " << filename
            << exit(FatalIOError);
    }

    OFstream& os = osPtr();
    bool ok = io.writeHeader(os, featureEdgeMesh::typeName);

    if (!ok)
    {
        FatalIOErrorInFunction
        (
            os
        )   << "Cannot write header"
            << exit(FatalIOError);
    }

    write(os, mesh.points(), mesh.edges());

    // Check state of Ostream
    os.check("edgeMeshFormat::write(Ostream&)");
}


// ************************************************************************* //
