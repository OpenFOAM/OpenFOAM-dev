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

Usage
    \b surfaceTransformPoints "\<transformations\>" \<input\> \<output\>
    Supported transformations:
      - \par translate=<translation vector>
        Translational transformation by given vector
      - \par rotate=(\<n1 vector\> \<n2 vector\>)
        Rotational transformation from unit vector n1 to n2
      - \par Rx=\<angle [deg] about x-axis\>
        Rotational transformation by given angle about x-axis
      - \par Ry=\<angle [deg] about y-axis\>
        Rotational transformation by given angle about y-axis
      - \par Rz=\<angle [deg] about z-axis\>
        Rotational transformation by given angle about z-axis
      - \par Ra=\<axis vector\> \<angle [deg] about axis\>
        Rotational transformation by given angle about given axis
      - \par scale=\<x-y-z scaling vector\>
        Anisotropic scaling by the given vector in the x, y, z
        coordinate directions

    Example usage:
        surfaceTransformPoints \
            "translate=(-0.586 0 -0.156), \
            Ry=3.485, \
            translate=(0.586 0 0.156)" \
            constant/geometry/w3_orig.stl constant/geometry/w3.stl

See also
    Foam::transformer
    transformPoints

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "unitConversion.H"
#include "MeshedSurfaces.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "removeCaseOptions.H"

    const wordList supportedTransformations
    (
        {"translate", "rotate", "Rx", "Ry", "Rz", "Ra", "scale"}
    );

    {
        OStringStream supportedTransformationsStr;
        supportedTransformationsStr << supportedTransformations << endl;

        argList::addNote
        (
            "Transforms a surface e.g.\n"
            "surfaceTransformPoints "
            "\"translate=(-0.586 0 -0.156), "
            "Ry=3.485, "
            "translate=(0.586 0 0.156)\" "
            "surf.stl tranformedSurf.obj\n\n"
            "Supported transformations " + supportedTransformationsStr.str()
        );
    }

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

    transformer transforms;

    forAll(transformations, i)
    {
        if (transformations[i].first() == "translate")
        {
            const vector v(IStringStream(transformations[i].second())());
            transforms = transformer::translation(v) & transforms;
        }
        else if (transformations[i].first() == "rotate")
        {
            Pair<vector> n1n2(IStringStream(transformations[i].second())());

            n1n2[0] /= mag(n1n2[0]);
            n1n2[1] /= mag(n1n2[1]);

            transforms =
                transformer::rotation(rotationTensor(n1n2[0], n1n2[1]))
              & transforms;
        }
        else if (transformations[i].first() == "Rx")
        {
            const scalar a
            (
                readScalar(IStringStream(transformations[i].second())())
            );
            transforms = transformer::rotation(Rx(degToRad(a))) & transforms;
        }
        else if (transformations[i].first() == "Ry")
        {
            const scalar a
            (
                readScalar(IStringStream(transformations[i].second())())
            );
            transforms = transformer::rotation(Ry(degToRad(a))) & transforms;
        }
        else if (transformations[i].first() == "Rz")
        {
            const scalar a
            (
                readScalar(IStringStream(transformations[i].second())())
            );
            transforms = transformer::rotation(Rz(degToRad(a))) & transforms;
        }
        else if (transformations[i].first() == "Ra")
        {
            IStringStream istr(transformations[i].second());
            const vector v(istr);
            const scalar a(readScalar(istr));
            transforms = transformer::rotation(Ra(v, degToRad(a))) & transforms;
        }
        else if (transformations[i].first() == "scale")
        {
            const vector v(IStringStream(transformations[i].second())());
            transforms =
                transformer::scaling(diagTensor(v.x(), v.y(), v.z()))
              & transforms;
        }
        else
        {
            args.printUsage();
            FatalErrorInFunction
                << "Unknown transformation " << transformations[i].first()
                << exit(FatalError);
        }
    }

    meshedSurface surf1(surfFileName);
    pointField points(surf1.points());
    transforms.transformPosition(points, points);
    surf1.movePoints(points);
    surf1.write(outFileName);

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
