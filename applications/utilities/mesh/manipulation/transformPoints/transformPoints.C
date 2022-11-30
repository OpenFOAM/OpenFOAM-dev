/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2022 OpenFOAM Foundation
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
    transformPoints

Description
    Transform (translate, rotate, scale) the mesh points, and optionally also
    any vector and tensor fields.

Usage
    \b transformPoints "\<transformations\>" [OPTION]

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

    Options:
      - \par -rotateFields \n
        Additionally transform vector and tensor fields.
      - \par -pointSet \<name\> \n
        Only transform points in the given point set.

    Example usage:
        transformPoints \
            "translate=(-0.05 -0.05 0), \
            Rz=45, \
            translate=(0.05 0.05 0)"

See also
    Foam::transformer
    surfaceTransformPoints

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "fvMesh.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "ReadFields.H"
#include "pointFields.H"
#include "pointSet.H"
#include "transformField.H"
#include "transformGeometricField.H"
#include "unitConversion.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class GeoField>
void readAndRotateFields
(
    PtrList<GeoField>& flds,
    const fvMesh& mesh,
    const tensor& T,
    const IOobjectList& objects
)
{
    ReadFields(mesh, objects, flds);
    forAll(flds, i)
    {
        Info<< "Transforming " << flds[i].name() << endl;
        dimensionedTensor dimT("t", flds[i].dimensions(), T);
        transform(flds[i], dimT, flds[i]);
    }
}


void rotateFields(const argList& args, const Time& runTime, const tensor& T)
{
    #include "createNamedMesh.H"

    // Read objects in time directory
    IOobjectList objects(mesh, runTime.name());

    // Read vol fields.
    PtrList<volScalarField> vsFlds;
    readAndRotateFields(vsFlds, mesh, T, objects);
    PtrList<volVectorField> vvFlds;
    readAndRotateFields(vvFlds, mesh, T, objects);
    PtrList<volSphericalTensorField> vstFlds;
    readAndRotateFields(vstFlds, mesh, T, objects);
    PtrList<volSymmTensorField> vsymtFlds;
    readAndRotateFields(vsymtFlds, mesh, T, objects);
    PtrList<volTensorField> vtFlds;
    readAndRotateFields(vtFlds, mesh, T, objects);

    // Read surface fields.
    PtrList<surfaceScalarField> ssFlds;
    readAndRotateFields(ssFlds, mesh, T, objects);
    PtrList<surfaceVectorField> svFlds;
    readAndRotateFields(svFlds, mesh, T, objects);
    PtrList<surfaceSphericalTensorField> sstFlds;
    readAndRotateFields(sstFlds, mesh, T, objects);
    PtrList<surfaceSymmTensorField> ssymtFlds;
    readAndRotateFields(ssymtFlds, mesh, T, objects);
    PtrList<surfaceTensorField> stFlds;
    readAndRotateFields(stFlds, mesh, T, objects);

    mesh.write();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Transforms a mesh by translation, rotation and/or scaling.\n"
        "The <transformations> are listed comma-separated in a string "
        "and executed in sequence.\n\n"
        "transformations:\n"
        "  translate=<vector>        "
        "translation by vector, e.g. (1 2 3)\n"
        "  rotate=(<n1> <n2>)        "
        "rotation from unit vector n1 to n2\n"
        "  Rx=<angle>                "
        "rotation by given angle [deg], e.g. 90, about x-axis\n"
        "  Ry=<angle>                "
        "rotation by given angle [deg] about y-axis\n"
        "  Rz=<angle>                "
        "rotation by given angle [deg] about z-axis\n"
        "  Ra=<axis vector> <angle>  "
        "rotation by given angle [deg] about specified axis\n"
        "  scale=<vector>            "
        "scale by factors from vector in x, y, z directions,\n"
        "                            "
        "e.g. (0.001 0.001 0.001) to scale from mm to m\n\n"
        "example:\n"
        "  transformPoints "
        "\"translate=(1.2 0 0), Rx=90, translate=(-1.2 0 0)\""
    );

    argList::validArgs.append("transformations");

    argList::addBoolOption
    (
        "rotateFields",
        "transform vector and tensor fields"
    );

    argList::addOption
    (
        "pointSet",
        "pointSet",
        "Point set to limit the transformation to"
    );

    #include "addRegionOption.H"
    #include "addAllRegionsOption.H"
    #include "setRootCase.H"
    #include "createTime.H"

    const string transformationString(args[1]);

    #include "createTransforms.H"

    #include "setRegionNames.H"

    const bool doRotateFields = args.optionFound("rotateFields");

    word pointSetName = word::null;
    const bool doPointSet = args.optionReadIfPresent("pointSet", pointSetName);

    if (doRotateFields && doPointSet)
    {
        FatalErrorInFunction
            << "Rotation of fields across the entire mesh, and limiting the "
            << "transformation of points to a set, cannot be done "
            << "simultaneously" << exit(FatalError);
    }

    forAll(regionNames, regioni)
    {
        const word& regionName = regionNames[regioni];

        const word& regionDir =
            regionName == polyMesh::defaultRegion
          ? word::null
          : regionName;

        const fileName meshDir(regionDir/polyMesh::meshSubDir);

        pointIOField points
        (
            IOobject
            (
                "points",
                runTime.findInstance(meshDir, "points"),
                meshDir,
                runTime,
                IOobject::MUST_READ,
                IOobject::NO_WRITE,
                false
            )
        );

        if (doPointSet)
        {
            const labelList setPointIDs
            (
                pointSet
                (
                    IOobject
                    (
                        pointSetName,
                        runTime.findInstance(meshDir/"sets", word::null),
                        polyMesh::meshSubDir/"sets",
                        runTime,
                        IOobject::MUST_READ,
                        IOobject::NO_WRITE,
                        false
                    )
                ).toc()
            );

            pointField setPoints(UIndirectList<point>(points, setPointIDs));

            transforms.transformPosition(setPoints, setPoints);

            UIndirectList<point>(points, setPointIDs) = setPoints;
        }
        else
        {
            transforms.transformPosition(points, points);
        }

        if (doRotateFields)
        {
            rotateFields(args, runTime, transforms.T());
        }

        // Set the precision of the points data to 10
        IOstream::defaultPrecision(max(10u, IOstream::defaultPrecision()));

        Info<< "Writing points into directory " << points.path() << nl << endl;
        points.write();
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
