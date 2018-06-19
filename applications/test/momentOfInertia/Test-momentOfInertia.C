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
    momentOfInertiaTest

Description
    Calculates the inertia tensor and principal axes and moments of a
    test face, tetrahedron and cell.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "polyMesh.H"
#include "ListOps.H"
#include "face.H"
#include "tetPointRef.H"
#include "triFaceList.H"
#include "OFstream.H"
#include "meshTools.H"
#include "momentOfInertia.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

using namespace Foam;

int main(int argc, char *argv[])
{
    argList::addOption
    (
        "cell",
        "label",
        "cell to use for inertia calculation, defaults to 0"
    );

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createPolyMesh.H"

    scalar density = 1.0;

    {
        label nPts = 6;

        pointField pts(nPts);

        pts[0] = point(4.495, 3.717, -4.112);
        pts[1] = point(4.421, 3.932, -4.112);
        pts[2] = point(4.379, 4.053, -4.112);
        pts[3] = point(4.301, 4.026, -4.300);
        pts[4] = point(4.294, 4.024, -4.317);
        pts[5] = point(4.409, 3.687, -4.317);

        face f(identity(nPts));

        point Cf = f.centre(pts);

        tensor J = Zero;

        J = f.inertia(pts, Cf, density);

        vector eVal = eigenValues(J);

        tensor eVec = eigenVectors(J);

        Info<< nl << "Inertia tensor of test face " << J << nl
            << "eigenValues (principal moments) " << eVal << nl
            << "eigenVectors (principal axes) " << eVec
            << endl;

        OFstream str("momentOfInertiaTestFace.obj");

        Info<< nl << "Writing test face and scaled principal axes to "
            << str.name() << endl;

        forAll(pts, ptI)
        {
            meshTools::writeOBJ(str, pts[ptI]);
        }

        str << "l";

        forAll(f, fI)
        {
            str << ' ' << fI + 1;
        }

        str << " 1" << endl;

        scalar scale = mag(Cf - pts[f[0]])/eVal.component(findMin(eVal));

        meshTools::writeOBJ(str, Cf);
        meshTools::writeOBJ(str, Cf + scale*eVal.x()*eVec.x());
        meshTools::writeOBJ(str, Cf + scale*eVal.y()*eVec.y());
        meshTools::writeOBJ(str, Cf + scale*eVal.z()*eVec.z());

        for (label i = nPts + 1; i < nPts + 4; i++)
        {
            str << "l " << nPts + 1 << ' ' << i + 1 << endl;
        }
    }

    {
        label nPts = 4;

        pointField pts(nPts);

        pts[0] = point(0, 0, 0);
        pts[1] = point(1, 0, 0);
        pts[2] = point(0.5, 1, 0);
        pts[3] = point(0.5, 0.5, 1);

        tetPointRef tet(pts[0], pts[1], pts[2], pts[3]);

        triFaceList tetFaces(4);

        tetFaces[0] = triFace(0, 2, 1);
        tetFaces[1] = triFace(1, 2, 3);
        tetFaces[2] = triFace(0, 3, 2);
        tetFaces[3] = triFace(0, 1, 3);

        scalar m = 0.0;
        vector cM = Zero;
        tensor J = Zero;

        momentOfInertia::massPropertiesSolid(pts, tetFaces, density, m, cM, J);

        vector eVal = eigenValues(J);

        tensor eVec = eigenVectors(J);

        Info<< nl
            << "Mass of tetrahedron " << m << nl
            << "Centre of mass of tetrahedron " << cM << nl
            << "Inertia tensor of tetrahedron " << J << nl
            << "eigenValues (principal moments) " << eVal << nl
            << "eigenVectors (principal axes) " << eVec
            << endl;

        OFstream str("momentOfInertiaTestTet.obj");

        Info<< nl << "Writing test tetrahedron and scaled principal axes to "
            << str.name() << endl;

        forAll(pts, ptI)
        {
            meshTools::writeOBJ(str, pts[ptI]);
        }

        forAll(tetFaces, tFI)
        {
            const triFace& f = tetFaces[tFI];

            str << "l";

            forAll(f, fI)
            {
                str << ' ' << f[fI] + 1;
            }

            str << ' ' << f[0] + 1 << endl;
        }

        scalar scale = mag(cM - pts[0])/eVal.component(findMin(eVal));

        meshTools::writeOBJ(str, cM);
        meshTools::writeOBJ(str, cM + scale*eVal.x()*eVec.x());
        meshTools::writeOBJ(str, cM + scale*eVal.y()*eVec.y());
        meshTools::writeOBJ(str, cM + scale*eVal.z()*eVec.z());

        for (label i = nPts + 1; i < nPts + 4; i++)
        {
            str << "l " << nPts + 1 << ' ' << i + 1 << endl;
        }
    }

    {
        const label celli = args.optionLookupOrDefault("cell", 0);

        tensorField mI(momentOfInertia::meshInertia(mesh));

        tensor& J = mI[celli];

        vector eVal = eigenValues(J);

        Info<< nl
            << "Inertia tensor of cell " << celli << " " << J << nl
            << "eigenValues (principal moments) " << eVal << endl;

        J /= cmptMax(eVal);

        tensor eVec = eigenVectors(J);

        Info<< "eigenVectors (principal axes, from normalised inertia) " << eVec
            << endl;

        OFstream str("cell_" + name(celli) + "_inertia.obj");

        Info<< nl << "Writing scaled principal axes of cell " << celli << " to "
            << str.name() << endl;

        const point& cC = mesh.cellCentres()[celli];

        scalar scale = mag
        (
            (cC - mesh.faceCentres()[mesh.cells()[celli][0]])
           /eVal.component(findMin(eVal))
        );

        meshTools::writeOBJ(str, cC);
        meshTools::writeOBJ(str, cC + scale*eVal.x()*eVec.x());
        meshTools::writeOBJ(str, cC + scale*eVal.y()*eVec.y());
        meshTools::writeOBJ(str, cC + scale*eVal.z()*eVec.z());

        for (label i = 1; i < 4; i++)
        {
            str << "l " << 1 << ' ' << i + 1 << endl;
        }
    }

    Info<< nl << "End" << nl << endl;

    return 0;
}


// ************************************************************************* //
