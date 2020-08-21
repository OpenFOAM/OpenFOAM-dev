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

#include <fstream>
#include <iostream>

using std::ofstream;
using std::ios;

#include "Time.H"
#include "fluentFvMesh.H"
#include "primitiveMesh.H"
#include "wallFvPatch.H"
#include "symmetryPlaneFvPatch.H"
#include "symmetryFvPatch.H"
#include "cellModeller.H"
#include "OSspecific.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fluentFvMesh::fluentFvMesh(const IOobject& io)
:
    fvMesh(io)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fluentFvMesh::writeFluentMesh() const
{
    // make a directory called proInterface in the case
    mkDir(time().rootPath()/time().caseName()/"fluentInterface");

    // open a file for the mesh
    ofstream fluentMeshFile
    (
        (
            time().rootPath()/
            time().caseName()/
            "fluentInterface"/
            time().caseName() + ".msh"
        ).c_str()
    );

    Info<< "Writing Header" << endl;

    fluentMeshFile
        << "(0 \"FOAM to Fluent Mesh File\")" << std::endl << std::endl
        << "(0 \"Dimension:\")" << std::endl
        << "(2 3)" << std::endl << std::endl
        << "(0 \"Grid dimensions:\")" << std::endl;

    // Writing number of points
    fluentMeshFile
            << "(10 (0 1 ";

    // Writing hex
    fluentMeshFile.setf(ios::hex, ios::basefield);

    fluentMeshFile
        << nPoints() << " 0 3))" << std::endl;

    // Writing number of cells
    fluentMeshFile
        << "(12 (0 1 "
        << nCells() << " 0 0))" << std::endl;

    // Writing number of faces
    label nFcs = nFaces();

    fluentMeshFile
            << "(13 (0 1 ";

    // Still writing hex
    fluentMeshFile
        << nFcs << " 0 0))" << std::endl << std::endl;

    // Return to dec
    fluentMeshFile.setf(ios::dec, ios::basefield);

    // Writing points
    fluentMeshFile
            << "(10 (1 1 ";

    fluentMeshFile.setf(ios::hex, ios::basefield);
    fluentMeshFile
        << nPoints() << " 1 3)"
        << std::endl << "(" << std::endl;

    fluentMeshFile.precision(10);
    fluentMeshFile.setf(ios::scientific);

    const pointField& p = points();

    forAll(p, pointi)
    {
        fluentMeshFile
            << "    "
            << p[pointi].x() << " "
            << p[pointi].y()
            << " " << p[pointi].z() << std::endl;
    }

    fluentMeshFile
        << "))" << std::endl << std::endl;

    const labelUList& own = owner();
    const labelUList& nei = neighbour();

    const faceList& fcs = faces();

    // Writing (mixed) internal faces
    fluentMeshFile
        << "(13 (2 1 "
        << own.size() << " 2 0)" << std::endl << "(" << std::endl;

    forAll(own, facei)
    {
        const labelList& l = fcs[facei];

        fluentMeshFile << "    ";

        fluentMeshFile << l.size() << " ";

        forAll(l, lI)
        {
            fluentMeshFile << l[lI] + 1 << " ";
        }

        fluentMeshFile << nei[facei] + 1 << " ";
        fluentMeshFile << own[facei] + 1 << std::endl;
    }

    fluentMeshFile << "))" << std::endl;

    label nWrittenFaces = own.size();

    // Writing boundary faces
    forAll(boundary(), patchi)
    {
        const faceUList& patchFaces = boundaryMesh()[patchi];

        const labelList& patchFaceCells =
            boundaryMesh()[patchi].faceCells();

        // The face group will be offset by 10 from the patch label

        // Write header
        fluentMeshFile
            << "(13 (" << patchi + 10 << " " << nWrittenFaces + 1
            << " " << nWrittenFaces + patchFaces.size() << " ";

        nWrittenFaces += patchFaces.size();

        // Write patch type
        if (isA<wallFvPatch>(boundary()[patchi]))
        {
            fluentMeshFile << 3;
        }
        else if
        (
            isA<symmetryPlaneFvPatch>(boundary()[patchi])
         || isA<symmetryFvPatch>(boundary()[patchi])
        )
        {
            fluentMeshFile << 7;
        }
        else
        {
            fluentMeshFile << 4;
        }

        fluentMeshFile
            <<" 0)" << std::endl << "(" << std::endl;

        forAll(patchFaces, facei)
        {
            const labelList& l = patchFaces[facei];

            fluentMeshFile << "    ";

            fluentMeshFile << l.size() << " ";

            // Note: In Fluent, all boundary faces point inwards, which is
            // opposite from the OpenFOAM convention.
            // Turn them around on printout
            forAllReverse (l, lI)
            {
                fluentMeshFile << l[lI] + 1 << " ";
            }

            fluentMeshFile << patchFaceCells[facei] + 1 << " 0" << std::endl;
        }

        fluentMeshFile << "))" << std::endl;
    }

    // Writing cells
    fluentMeshFile
        << "(12 (1 1 "
        << nCells() << " 1 0)(" << std::endl;

    const cellModel& hex = *(cellModeller::lookup("hex"));
    const cellModel& prism = *(cellModeller::lookup("prism"));
    const cellModel& pyr = *(cellModeller::lookup("pyr"));
    const cellModel& tet = *(cellModeller::lookup("tet"));

    const cellShapeList& cells = cellShapes();

    bool hasWarned = false;

    forAll(cells, celli)
    {
        if (cells[celli].model() == tet)
        {
            fluentMeshFile << " " << 2;
        }
        else if (cells[celli].model() == hex)
        {
            fluentMeshFile << " " << 4;
        }
        else if (cells[celli].model() == pyr)
        {
            fluentMeshFile << " " << 5;
        }
        else if (cells[celli].model() == prism)
        {
            fluentMeshFile << " " << 6;
        }
        else
        {
            if (!hasWarned)
            {
                hasWarned = true;

                WarningInFunction
                    << "foamMeshToFluent: cell shape for cell "
                    << celli << " only supported by Fluent polyhedral meshes."
                    << nl
                    << "    Suppressing any further messages for polyhedral"
                    << " cells." << endl;
            }
            fluentMeshFile << " " << 7;
        }
    }

    fluentMeshFile << ")())" << std::endl;

    // Return to dec
    fluentMeshFile.setf(ios::dec, ios::basefield);

    // Writing patch types
    fluentMeshFile << "(39 (1 fluid fluid-1)())" << std::endl;
    fluentMeshFile << "(39 (2 interior interior-1)())" << std::endl;

    // Writing boundary patch types
    forAll(boundary(), patchi)
    {
        fluentMeshFile
            << "(39 (" << patchi + 10 << " ";

        // Write patch type
        if (isA<wallFvPatch>(boundary()[patchi]))
        {
            fluentMeshFile << "wall ";
        }
        else if
        (
            isA<symmetryPlaneFvPatch>(boundary()[patchi])
         || isA<symmetryFvPatch>(boundary()[patchi])
        )
        {
            fluentMeshFile << "symmetry ";
        }
        else
        {
            fluentMeshFile << "pressure-outlet ";
        }

        fluentMeshFile
            << boundary()[patchi].name() << ")())" << std::endl;
    }
}


// ************************************************************************* //
