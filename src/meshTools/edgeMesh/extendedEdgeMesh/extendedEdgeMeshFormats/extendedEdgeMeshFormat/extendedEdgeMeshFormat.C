/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
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

#include "extendedEdgeMeshFormat.H"
#include "IFstream.H"
#include "Time.H"
#include "extendedFeatureEdgeMesh.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fileFormats::extendedEdgeMeshFormat::extendedEdgeMeshFormat
(
    const fileName& filename
)
{
    read(filename);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::fileFormats::extendedEdgeMeshFormat::read
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

    if (!io.typeHeaderOk<extendedFeatureEdgeMesh>(false))
    {
        FatalErrorInFunction
            << "Cannot read file " << filename
            << exit(FatalError);
    }

    const fileName fName(typeFilePath<extendedFeatureEdgeMesh>(io));

    autoPtr<IFstream> isPtr(new IFstream(fName));
    bool ok = false;
    if (isPtr().good())
    {
        Istream& is = isPtr();
        ok = io.readHeader(is);

        if (ok)
        {
            // Use extendedEdgeMesh IO
            is >> *this;
            ok = is.good();
        }
    }

    return ok;
}


// ************************************************************************* //
