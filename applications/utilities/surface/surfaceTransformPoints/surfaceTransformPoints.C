/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2021 OpenFOAM Foundation
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

Application
    surfaceTransformPoints

Description
    Transform (translate, rotate, scale) a surface.

    The rollPitchYaw option takes three angles (degrees):
    - roll (rotation about x) followed by
    - pitch (rotation about y) followed by
    - yaw (rotation about z)

    The yawPitchRoll does yaw followed by pitch followed by roll.

Usage
    \b surfaceTransformPoints "\<transformations\>" \<input\> \<output\>

    Example usage:
        surfaceTransformPoints \
            "translate=(-0.586 0 -0.156), \
            rollPitchYaw=(0 -3.485 0), \
            translate=(0.586 0 0.156)" \
            constant/geometry/w3_orig.stl constant/geometry/w3.stl

See also
    transformPoints

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "quaternion.H"
#include "unitConversion.H"
#include "MeshedSurfaces.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "removeCaseOptions.H"

    argList::addNote
    (
        "Transform (translate, rotate, scale) a surface."
    );
    argList::noParallel();
    argList::validArgs.append("transformations");
    argList::validArgs.append("surface file");
    argList::validArgs.append("output surface file");

    argList args(argc, argv);

    const string transformationString(args[1]);
    const fileName surfFileName(args[2]);
    const fileName outFileName(args[3]);

    Info<< "Reading surf from " << surfFileName << " ..." << nl
        << "Writing surf to " << outFileName << " ..." << endl;

    wordReList simpleTransformations;
    List<Tuple2<word, string>> transformations;
    dictArgList(transformationString, simpleTransformations, transformations);

    meshedSurface surf1(surfFileName);

    pointField points(surf1.points());

    forAll(transformations, i)
    {
        if (transformations[i].first() == "translate")
        {
            const vector v(IStringStream(transformations[i].second())());
            Info<< "Translating points by " << v << endl;
            points += v;
        }
        else if (transformations[i].first() == "rotate")
        {
            Pair<vector> n1n2(IStringStream(transformations[i].second())());

            n1n2[0] /= mag(n1n2[0]);
            n1n2[1] /= mag(n1n2[1]);

            const tensor T = rotationTensor(n1n2[0], n1n2[1]);

            Info<< "Rotating points by " << T << endl;

            points = transform(T, points);
        }
        else if (transformations[i].first() == "rollPitchYaw")
        {
            const vector v(IStringStream(transformations[i].second())());

            Info<< "Rotating points by "
                << " roll = " << v.x()
                << ", pitch = " << v.y()
                << ", yaw = " << v.z() << nl;

            const quaternion R(quaternion::rotationSequence::XYZ, degToRad(v));
            points = transform(R, points);
        }
        else if (transformations[i].first() == "yawPitchRoll")
        {
            const vector v(IStringStream(transformations[i].second())());

            Info<< "Rotating points by "
                << " yaw   " << v.x()
                << ", pitch " << v.y()
                << ", roll  " << v.z() << nl;

            const quaternion R
            (
                quaternion::rotationSequence::ZYX,
                degToRad(vector(v.z(), v.y(), v.x()))
            );
            points = transform(R, points);
        }
        else if (transformations[i].first() == "scale")
        {
            const vector v(IStringStream(transformations[i].second())());

            Info<< "Scaling points by " << v << endl;

            points.replace(vector::X, v.x()*points.component(vector::X));
            points.replace(vector::Y, v.y()*points.component(vector::Y));
            points.replace(vector::Z, v.z()*points.component(vector::Z));
        }
    }

    surf1.movePoints(points);
    surf1.write(outFileName);

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
