/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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

#include "surfaceMeshWriter.H"
#include "writeFuns.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::surfaceMeshWriter::surfaceMeshWriter
(
    const bool binary,
    const indirectPrimitivePatch& pp,
    const word& name,
    const fileName& fName
)
:
    binary_(binary),
    pp_(pp),
    fName_(fName),
    os_(fName.c_str())
{
    // Write header
    writeFuns::writeHeader(os_, binary_, name);

    os_ << "DATASET POLYDATA" << std::endl;

    // Write topology
    label nFaceVerts = 0;

    forAll(pp, faceI)
    {
        nFaceVerts += pp[faceI].size() + 1;
    }

    os_ << "POINTS " << pp.nPoints() << " float" << std::endl;

    DynamicList<floatScalar> ptField(3*pp.nPoints());
    writeFuns::insert(pp.localPoints(), ptField);
    writeFuns::write(os_, binary, ptField);


    os_ << "POLYGONS " << pp.size() << ' ' << nFaceVerts << std::endl;

    DynamicList<label> vertLabels(nFaceVerts);

    forAll(pp, faceI)
    {
        const face& f = pp.localFaces()[faceI];

        vertLabels.append(f.size());
        writeFuns::insert(f, vertLabels);
    }
    writeFuns::write(os_, binary_, vertLabels);
}


// ************************************************************************* //
