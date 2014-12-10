/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

#include "dxSurfaceWriter.H"

#include "OFstream.H"
#include "OSspecific.H"

#include "makeSurfaceWriterMethods.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    makeSurfaceWriterType(dxSurfaceWriter);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::dxSurfaceWriter::writeGeometry
(
    Ostream& os,
    const pointField& points,
    const faceList& faces
)
{
    // Write vertex coordinates

    os  << "# The irregular positions" << nl
        << "object 1 class array type float rank 1 shape 3 items "
        << points.size() << " data follows" << nl;

    forAll(points, pointI)
    {
        const point& pt = points[pointI];

        os  << float(pt.x()) << ' ' << float(pt.y()) << ' ' << float(pt.z())
            << nl;
    }
    os  << nl;

    // Write triangles
    os  << "# The irregular connections (triangles)" << nl
        << "object 2 class array type int rank 1 shape 3 items "
        << faces.size() << " data follows" << nl;

    forAll(faces, faceI)
    {
        const face& f = faces[faceI];

        if (f.size() != 3)
        {
            FatalErrorIn
            (
                "writeGeometry(Ostream&, const pointField&, const faceList&)"
            )   << "Face " << faceI << " vertices " << f
                << " is not a triangle."
                << exit(FatalError);
        }

        os << f[0] << ' ' << f[1] << ' ' << f[2] << nl;
    }
    os << "attribute \"element type\" string \"triangles\"" << nl
       << "attribute \"ref\" string \"positions\"" << nl << nl;
}


void Foam::dxSurfaceWriter::writeTrailer(Ostream& os, const bool isNodeValues)
{
    if (isNodeValues)
    {
        os  << nl << "attribute \"dep\" string \"positions\""
            << nl << nl;
    }
    else
    {
        os  << nl << "attribute \"dep\" string \"connections\""
            << nl << nl;
    }

    os  << "# the field, with three components: \"positions\","
        << " \"connections\", and \"data\"" << nl
        << "object \"irregular positions irregular "
        << "connections\" class field"
        << nl
        << "component \"positions\" value 1" << nl
        << "component \"connections\" value 2" << nl
        << "component \"data\" value 3" << nl;

    os  << "end" << endl;
}


namespace Foam
{
    template<>
    void Foam::dxSurfaceWriter::writeData
    (
        Ostream& os,
        const Field<scalar>& values
    )
    {
        os  << "object 3 class array type float rank 0 items "
            << values.size() << " data follows" << nl;

        forAll(values, elemI)
        {
            os << float(values[elemI]) << nl;
        }
    }


    template<>
    void Foam::dxSurfaceWriter::writeData
    (
        Ostream& os,
        const Field<vector>& values
    )
    {
        os  << "object 3 class array type float rank 1 shape 3 items "
            << values.size() << " data follows" << nl;

        forAll(values, elemI)
        {
            os  << float(values[elemI].x()) << ' '
                << float(values[elemI].y()) << ' '
                << float(values[elemI].z()) << nl;
        }
    }


    template<>
    void Foam::dxSurfaceWriter::writeData
    (
        Ostream& os,
        const Field<sphericalTensor>& values
    )
    {
        os  << "object 3 class array type float rank 0 items "
            << values.size() << " data follows" << nl;

        forAll(values, elemI)
        {
            os << float(values[elemI][0]) << nl;
        }
    }


    template<>
    void Foam::dxSurfaceWriter::writeData
    (
        Ostream& os,
        const Field<symmTensor>& values
    )
    {
        os  << "object 3 class array type float rank 2 shape 3 items "
            << values.size() << " data follows" << nl;

        forAll(values, elemI)
        {
            const symmTensor& t = values[elemI];

            os  << float(t.xx()) << ' ' << float(t.xy()) << ' ' << float(t.xz())
                << float(t.xy()) << ' ' << float(t.yy()) << ' ' << float(t.yz())
                << float(t.xz()) << ' ' << float(t.yz()) << ' ' << float(t.zz())
                << nl;
        }
    }


    // Write Field<tensor> in DX format
    template<>
    inline void Foam::dxSurfaceWriter::writeData
    (
        Ostream& os,
        const Field<tensor>& values
    )
    {
        os  << "object 3 class array type float rank 2 shape 3 items "
            << values.size() << " data follows" << nl;

        forAll(values, elemI)
        {
            const tensor& t = values[elemI];

            os  << float(t.xx()) << ' ' << float(t.xy()) << ' ' << float(t.xz())
                << float(t.yx()) << ' ' << float(t.yy()) << ' ' << float(t.yz())
                << float(t.zx()) << ' ' << float(t.zy()) << ' ' << float(t.zz())
                << nl;
        }
    }
}


// arbitrary field
template<class Type>
inline void Foam::dxSurfaceWriter::writeData
(
    Ostream& os,
    const Field<Type>& values
)
{
    os  << "object 3 class array type float rank 0 items "
        << values.size() << " data follows" << nl;

    forAll(values, elemI)
    {
        os << float(0.0) << nl;
    }
}


template<class Type>
void Foam::dxSurfaceWriter::writeTemplate
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

    OFstream os
    (
        outputDir/fieldName + '_' + surfaceName + ".dx"
    );

    if (verbose)
    {
        Info<< "Writing field " << fieldName << " to " << os.name() << endl;
    }

    writeGeometry(os, points, faces);
    writeData(os, values);
    writeTrailer(os, isNodeValues);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::dxSurfaceWriter::dxSurfaceWriter()
:
    surfaceWriter()
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::dxSurfaceWriter::~dxSurfaceWriter()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// create write methods
defineSurfaceWriterWriteFields(Foam::dxSurfaceWriter);


// ************************************************************************* //
