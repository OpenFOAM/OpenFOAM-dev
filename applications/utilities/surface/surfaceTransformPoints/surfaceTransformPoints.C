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

Application
    surfaceTransformPoints

Description
    Transform (scale/rotate) a surface.
    Like transformPoints but for surfaces.

    The rollPitchYaw option takes three angles (degrees):
    - roll (rotation about x) followed by
    - pitch (rotation about y) followed by
    - yaw (rotation about z)

    The yawPitchRoll does yaw followed by pitch followed by roll.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "OFstream.H"
#include "IFstream.H"
#include "boundBox.H"
#include "transformField.H"
#include "Pair.H"
#include "quaternion.H"
#include "mathematicalConstants.H"

#include "MeshedSurfaces.H"

using namespace Foam;
using namespace Foam::constant::mathematical;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "removeCaseOptions.H"

    argList::addNote
    (
        "Transform (scale/rotate) a surface. "
        "Like transformPoints but for surfaces."
    );
    argList::noParallel();
    argList::validArgs.append("surface file");
    argList::validArgs.append("output surface file");
    argList::addOption
    (
        "translate",
        "vector",
        "translate by the specified <vector> - eg, '(1 0 0)'"
    );
    argList::addOption
    (
        "rotate",
        "(vectorA vectorB)",
        "transform in terms of a rotation between <vectorA> and <vectorB> "
        "- eg, '( (1 0 0) (0 0 1) )'"
    );
    argList::addOption
    (
        "scale",
        "vector",
        "scale by the specified amount - eg, '(0.001 0.001 0.001)' for a "
        "uniform [mm] to [m] scaling"
    );
    argList::addOption
    (
        "rollPitchYaw",
        "vector",
        "transform in terms of '( roll pitch yaw )' in degrees"
    );
    argList::addOption
    (
        "yawPitchRoll",
        "vector",
        "transform in terms of '( yaw pitch roll )' in degrees"
    );
    argList args(argc, argv);

    const fileName surfFileName = args[1];
    const fileName outFileName  = args[2];

    Info<< "Reading surf from " << surfFileName << " ..." << nl
        << "Writing surf to " << outFileName << " ..." << endl;

    if (args.options().empty())
    {
        FatalErrorInFunction
            << "No options supplied, please use one or more of "
               "-translate, -rotate or -scale options."
            << exit(FatalError);
    }

    meshedSurface surf1(surfFileName);

    pointField points(surf1.points());

    vector v;
    if (args.optionReadIfPresent("translate", v))
    {
        Info<< "Translating points by " << v << endl;

        points += v;
    }

    if (args.optionFound("rotate"))
    {
        Pair<vector> n1n2
        (
            args.optionLookup("rotate")()
        );
        n1n2[0] /= mag(n1n2[0]);
        n1n2[1] /= mag(n1n2[1]);

        tensor T = rotationTensor(n1n2[0], n1n2[1]);

        Info<< "Rotating points by " << T << endl;

        points = transform(T, points);
    }
    else if (args.optionReadIfPresent("rollPitchYaw", v))
    {
        Info<< "Rotating points by" << nl
            << "    roll  " << v.x() << nl
            << "    pitch " << v.y() << nl
            << "    yaw   " << v.z() << nl;

        // Convert to radians
        v *= pi/180.0;

        quaternion R(quaternion::rotationSequence::XYZ, v);

        Info<< "Rotating points by quaternion " << R << endl;
        points = transform(R, points);
    }
    else if (args.optionReadIfPresent("yawPitchRoll", v))
    {
        Info<< "Rotating points by" << nl
            << "    yaw   " << v.x() << nl
            << "    pitch " << v.y() << nl
            << "    roll  " << v.z() << nl;


        // Convert to radians
        v *= pi/180.0;

        scalar yaw = v.x();
        scalar pitch = v.y();
        scalar roll = v.z();

        quaternion R = quaternion(vector(0, 0, 1), yaw);
        R *= quaternion(vector(0, 1, 0), pitch);
        R *= quaternion(vector(1, 0, 0), roll);

        Info<< "Rotating points by quaternion " << R << endl;
        points = transform(R, points);
    }

    if (args.optionReadIfPresent("scale", v))
    {
        Info<< "Scaling points by " << v << endl;

        points.replace(vector::X, v.x()*points.component(vector::X));
        points.replace(vector::Y, v.y()*points.component(vector::Y));
        points.replace(vector::Z, v.z()*points.component(vector::Z));
    }

    surf1.movePoints(points);
    surf1.write(outFileName);

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
