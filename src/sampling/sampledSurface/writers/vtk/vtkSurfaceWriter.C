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

#include "vtkSurfaceWriter.H"
#include "OFstream.H"
#include "OSspecific.H"
#include "makeSurfaceWriterMethods.H"
#include "vtkWriteOps.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    makeSurfaceWriterType(vtkSurfaceWriter);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::vtkSurfaceWriter::writeGeometry
(
    std::ostream& os,
    const pointField& points,
    const faceList& faces
) const
{
    const bool binary = (writeFormat_ == IOstream::BINARY);

    // VTK header
    vtkWriteOps::writeHeader(os, binary, "sampleSurface");
    os << "DATASET POLYDATA" << nl;

    // Write vertex coords
    os  << "POINTS " << points.size() << " float" << nl;

    List<floatScalar> po(points.size()*3);
    label ind = 0;
    forAll(points, pointi)
    {
        const point& pt = points[pointi];
        forAll(pt, cmpt)
        {
            po[ind++] = float(pt[cmpt]);
        }
    }
    vtkWriteOps::write(os, binary, po);

    // Write faces
    label nNodes = 0;
    forAll(faces, facei)
    {
        nNodes += faces[facei].size();
    }

    os  << "POLYGONS " << faces.size() << ' '
        << faces.size() + nNodes << nl;

    labelList polygons(faces.size() + nNodes);
    ind = 0;
    forAll(faces, facei)
    {
        const face& f = faces[facei];
        polygons[ind++] = f.size();
        forAll(f, fp)
        {
            polygons[ind++] = f[fp];
        }
    }
    vtkWriteOps::write(os, binary, polygons);
}


template<class Type>
void Foam::vtkSurfaceWriter::Write
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
    const bool binary = (writeFormat_ == IOstream::BINARY);

    if (!isDir(outputDir))
    {
        mkDir(outputDir);
    }

    const word filePath = outputDir/fieldName + '_' + surfaceName + ".vtk";

    ofstream os(filePath, std::ios::binary);

    if (debug)
    {
        Info<< "Writing field " << fieldName << " to " << filePath << endl;
    }

    writeGeometry(os, points, faces);

    // Write data
    if (isNodeValues)
    {
        os  << "POINT_DATA ";
    }
    else
    {
        os  << "CELL_DATA ";
    }

    os  << values.size() << nl
        << "FIELD attributes 1" << nl
        << fieldName << " ";

    const label nComp = pTraits<Type>::nComponents;

    os  << nComp << " " << values.size() << " float" << nl;

    List<floatScalar> vals(values.size()*nComp);
    label ind = 0;
    forAll(values, elemI)
    {
        for (direction cmpt=0; cmpt < nComp; ++cmpt)
        {
            vals[ind++] = component(values[elemI], cmpt);
        }
    }

    vtkWriteOps::write(os, binary, vals);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::vtkSurfaceWriter::vtkSurfaceWriter
(
    const IOstream::streamFormat writeFormat
)
:
    surfaceWriter(writeFormat)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::vtkSurfaceWriter::~vtkSurfaceWriter()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::vtkSurfaceWriter::write
(
    const fileName& outputDir,
    const fileName& surfaceName,
    const pointField& points,
    const faceList& faces
) const
{
    if (!isDir(outputDir))
    {
        mkDir(outputDir);
    }

    word filePath =  outputDir/surfaceName + ".vtk";
    ofstream os(filePath, std::ios::binary);

    if (debug)
    {
        Info<< "Writing geometry to " << filePath << endl;
    }

    writeGeometry(os, points, faces);
}


// Create write methods
defineSurfaceWriterWriteFields(Foam::vtkSurfaceWriter);


// ************************************************************************* //
