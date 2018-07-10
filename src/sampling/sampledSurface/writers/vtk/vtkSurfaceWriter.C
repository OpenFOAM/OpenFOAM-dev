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

#include "vtkSurfaceWriter.H"

#include "OFstream.H"
#include "OSspecific.H"

#include "makeSurfaceWriterMethods.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    makeSurfaceWriterType(vtkSurfaceWriter);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::vtkSurfaceWriter::writeGeometry
(
    Ostream& os,
    const pointField& points,
    const faceList& faces
)
{
    // header
    os
        << "# vtk DataFile Version 2.0" << nl
        << "sampleSurface" << nl
        << "ASCII" << nl
        << "DATASET POLYDATA" << nl;

    // Write vertex coords
    os  << "POINTS " << points.size() << " double" << nl;
    forAll(points, pointi)
    {
        const point& pt = points[pointi];
        os  << float(pt.x()) << ' '
            << float(pt.y()) << ' '
            << float(pt.z()) << nl;
    }
    os  << nl;


    // Write faces
    label nNodes = 0;
    forAll(faces, facei)
    {
        nNodes += faces[facei].size();
    }

    os  << "POLYGONS " << faces.size() << ' '
        << faces.size() + nNodes << nl;

    forAll(faces, facei)
    {
        const face& f = faces[facei];

        os  << f.size();
        forAll(f, fp)
        {
            os  << ' ' << f[fp];
        }
        os  << nl;
    }
}


namespace Foam
{

    template<>
    void Foam::vtkSurfaceWriter::writeData
    (
        Ostream& os,
        const Field<scalar>& values
    )
    {
        os  << "1 " << values.size() << " float" << nl;

        forAll(values, elemI)
        {
            if (elemI)
            {
                if (elemI % 10)
                {
                    os  << ' ';
                }
                else
                {
                    os  << nl;
                }
            }

            os  << float(values[elemI]);
        }
        os  << nl;
    }


    template<>
    void Foam::vtkSurfaceWriter::writeData
    (
        Ostream& os,
        const Field<vector>& values
    )
    {
        os  << "3 " << values.size() << " float" << nl;

        forAll(values, elemI)
        {
            const vector& v = values[elemI];
            os  << float(v[0]) << ' ' << float(v[1]) << ' ' << float(v[2])
                << nl;
        }
    }


    template<>
    void Foam::vtkSurfaceWriter::writeData
    (
        Ostream& os,
        const Field<sphericalTensor>& values
    )
    {
        os  << "1 " << values.size() << " float" << nl;

        forAll(values, elemI)
        {
            const sphericalTensor& v = values[elemI];
            os  << float(v[0]) << nl;
        }
    }


    template<>
    void Foam::vtkSurfaceWriter::writeData
    (
        Ostream& os,
        const Field<symmTensor>& values
    )
    {
        os  << "6 " << values.size() << " float" << nl;

        forAll(values, elemI)
        {
            const symmTensor& v = values[elemI];
            os  << float(v[0]) << ' ' << float(v[1]) << ' ' << float(v[2])
                << ' '
                << float(v[3]) << ' ' << float(v[4]) << ' ' << float(v[5])
                << nl;

        }
    }


    template<>
    void Foam::vtkSurfaceWriter::writeData
    (
        Ostream& os,
        const Field<tensor>& values
    )
    {
        os  << "9 " << values.size() << " float" << nl;

        forAll(values, elemI)
        {
            const tensor& v = values[elemI];
            os  << float(v[0]) << ' ' << float(v[1]) << ' ' << float(v[2])
                << ' '
                << float(v[3]) << ' ' << float(v[4]) << ' ' << float(v[5])
                << ' '
                << float(v[6]) << ' ' << float(v[7]) << ' ' << float(v[8])
                << nl;
        }
    }

}


// Write generic field in vtk format
template<class Type>
void Foam::vtkSurfaceWriter::writeData
(
    Ostream& os,
    const Field<Type>& values
)
{
    os  << "1 " << values.size() << " float" << nl;

    forAll(values, elemI)
    {
        os  << float(0) << nl;
    }
}


template<class Type>
void Foam::vtkSurfaceWriter::writeTemplate
(
    const fileName& outputDir,
    const fileName& surfaceName,
    const pointField& points,
    const faceList& faces,
    const word& fieldName,
    const Field<Type>& values,
    const bool isNodeValues,
    const bool verbose
) const
{
    if (!isDir(outputDir))
    {
        mkDir(outputDir);
    }

    OFstream os(outputDir/fieldName + '_' + surfaceName + ".vtk");

    if (verbose)
    {
        Info<< "Writing field " << fieldName << " to " << os.name() << endl;
    }

    writeGeometry(os, points, faces);

    // start writing data
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

    // Write data
    writeData(os, values);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::vtkSurfaceWriter::vtkSurfaceWriter()
:
    surfaceWriter()
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
    const faceList& faces,
    const bool verbose
) const
{
    if (!isDir(outputDir))
    {
        mkDir(outputDir);
    }

    OFstream os(outputDir/surfaceName + ".vtk");

    if (verbose)
    {
        Info<< "Writing geometry to " << os.name() << endl;
    }

    writeGeometry(os, points, faces);
}


// create write methods
defineSurfaceWriterWriteFields(Foam::vtkSurfaceWriter);


// ************************************************************************* //
