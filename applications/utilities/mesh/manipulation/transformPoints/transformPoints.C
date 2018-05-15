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

Application
    transformPoints

Description
    Transforms the mesh points in the polyMesh directory according to the
    translate, rotate and scale options.

Usage
    \b transformPoints [OPTION]

    Options:
      - \par -translate \<vector\> \n
        Translates the points by the given vector.

      - \par -rotate (\<vector\> \<vector\>) \n
        Rotates the points from the first vector to the second.

      - \par -yawPitchRoll (\<yawdegrees\> \<pitchdegrees\> \<rolldegrees\>) \n
        Alternative rotation specification:
            yaw (rotation about z)
            pitch (rotation about y)
            roll (rotation about x)

      - \par -rollPitchYaw (\<rolldegrees\> \<pitchdegrees\> \<yawdegrees\>) \n
        Alternative rotation specification:
            roll (rotation about x)
            pitch (rotation about y)
            yaw (rotation about z)

      - \par -rotateFields \n
        In combination with \a -rotate, \a -yawPitchRoll or \a -rollPitchYaw
        additionally transform vector and tensor fields.

      - \par -scale \<vector\> \n
        Scales the points by the given vector.

    Any or all of the three transformation option types may be specified and are
    processed in the above order.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "fvMesh.H"
#include "regionProperties.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "ReadFields.H"
#include "pointFields.H"
#include "transformField.H"
#include "transformGeometricField.H"
#include "mathematicalConstants.H"

using namespace Foam;
using namespace Foam::constant::mathematical;

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
    IOobjectList objects(mesh, runTime.timeName());

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
        "rollPitchYaw",
        "vector",
        "transform in terms of '(roll pitch yaw)' in degrees"
    );
    argList::addOption
    (
        "yawPitchRoll",
        "vector",
        "transform in terms of '(yaw pitch roll)' in degrees"
    );
    argList::addBoolOption
    (
        "rotateFields",
        "read and transform vector and tensor fields too"
    );
    argList::addOption
    (
        "scale",
        "vector",
        "scale by the specified amount - eg, '(0.001 0.001 0.001)' for a "
        "uniform [mm] to [m] scaling"
    );

    #include "addRegionOption.H"
    #include "addAllRegionsOption.H"
    #include "setRootCase.H"
    #include "createTime.H"

    const wordList regionNames(selectRegionNames(args, runTime));

    forAll(regionNames, regioni)
    {
        const word& regionName = regionNames[regioni];
        const word& regionDir = Foam::regionDir(regionName);

        fileName meshDir(regionDir/polyMesh::meshSubDir);

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

        const bool doRotateFields = args.optionFound("rotateFields");

        // this is not actually stringent enough:
        if (args.options().empty())
        {
            FatalErrorInFunction
                << "No options supplied, please use one or more of "
                   "-translate, -rotate or -scale options."
                << exit(FatalError);
        }

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

            if (doRotateFields)
            {
                rotateFields(args, runTime, T);
            }
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

            if (doRotateFields)
            {
                rotateFields(args, runTime, R.R());
            }
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

            if (doRotateFields)
            {
                rotateFields(args, runTime, R.R());
            }
        }

        if (args.optionReadIfPresent("scale", v))
        {
            Info<< "Scaling points by " << v << endl;

            points.replace(vector::X, v.x()*points.component(vector::X));
            points.replace(vector::Y, v.y()*points.component(vector::Y));
            points.replace(vector::Z, v.z()*points.component(vector::Z));
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
