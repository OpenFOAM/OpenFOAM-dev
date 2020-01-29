/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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

#include "foamSurfaceWriter.H"
#include "OFstream.H"
#include "OSspecific.H"
#include "makeSurfaceWriterMethods.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    makeSurfaceWriterType(foamSurfaceWriter);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
void Foam::foamSurfaceWriter::Write
(
    const fileName& outputDir,
    const fileName& surfaceName,
    const pointField& points,
    const faceList& faces,
    const word& fieldName,
    const Field<Type>& values,
    const bool isNodeValues
) const
{
    const fileName surfaceDir(outputDir/surfaceName);

    if (!isDir(surfaceDir))
    {
        mkDir(surfaceDir);
    }

    if (debug)
    {
        Info<< "Writing field " << fieldName << " to " << surfaceDir << endl;
    }

    // Geometry should already have been written
    // Values to separate directory (e.g. "scalarField/p")

    const fileName foamName(pTraits<Type>::typeName);
    const fileName valuesDir(surfaceDir/(foamName + Field<Type>::typeName));

    if (!isDir(valuesDir))
    {
        mkDir(valuesDir);
    }

    OFstream(valuesDir/fieldName, writeFormat_)()  << values;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::foamSurfaceWriter::foamSurfaceWriter
(
    const IOstream::streamFormat writeFormat
)
:
    surfaceWriter(writeFormat)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::foamSurfaceWriter::~foamSurfaceWriter()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::foamSurfaceWriter::write
(
    const fileName& outputDir,
    const fileName& surfaceName,
    const pointField& points,
    const faceList& faces
) const
{
    const fileName surfaceDir(outputDir/surfaceName);

    if (!isDir(surfaceDir))
    {
        mkDir(surfaceDir);
    }

    if (debug)
    {
        Info<< "Writing geometry to " << surfaceDir << endl;
    }

    // Points
    OFstream(surfaceDir/"points", writeFormat_)() << points;

    // Faces
    OFstream(surfaceDir/"faces", writeFormat_)() << faces;

    // Face centers. Not really necessary but very handy when reusing as inputs
    // for e.g. timeVaryingMapped bc.
    pointField faceCentres(faces.size(), Zero);

    forAll(faces, facei)
    {
        faceCentres[facei] = faces[facei].centre(points);
    }

    OFstream(surfaceDir/"faceCentres", writeFormat_)() << faceCentres;
}


// Create write methods
defineSurfaceWriterWriteFields(Foam::foamSurfaceWriter);


// ************************************************************************* //
