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

#include "patchWriter.H"
#include "vtkWriteFieldOps.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::patchWriter::patchWriter
(
    const vtkMesh& vMesh,
    const bool binary,
    const bool nearCellValue,
    const fileName& fName,
    const labelList& patchIDs
)
:
    vMesh_(vMesh),
    binary_(binary),
    nearCellValue_(nearCellValue),
    fName_(fName),
    patchIDs_(patchIDs),
    os_(fName.c_str())
{
    const fvMesh& mesh = vMesh_.mesh();
    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    // Write header
    if (patchIDs_.size() == 1)
    {
        vtkWriteOps::writeHeader(os_, binary_, patches[patchIDs_[0]].name());
    }
    else
    {
        vtkWriteOps::writeHeader(os_, binary_, "patches");
    }
    os_ << "DATASET POLYDATA" << std::endl;

    // Write topology
    nPoints_ = 0;
    nFaces_ = 0;
    label nFaceVerts = 0;

    forAll(patchIDs_, i)
    {
        const polyPatch& pp = patches[patchIDs_[i]];

        nPoints_ += pp.nPoints();
        nFaces_ += pp.size();

        forAll(pp, facei)
        {
            nFaceVerts += pp[facei].size() + 1;
        }
    }

    os_ << "POINTS " << nPoints_ << " float" << std::endl;

    DynamicList<floatScalar> ptField(3*nPoints_);

    forAll(patchIDs_, i)
    {
        const polyPatch& pp = patches[patchIDs_[i]];

        vtkWriteOps::insert(pp.localPoints(), ptField);
    }
    vtkWriteOps::write(os_, binary_, ptField);

    os_ << "POLYGONS " << nFaces_ << ' ' << nFaceVerts << std::endl;

    DynamicList<label> vertLabels(nFaceVerts);

    label offset = 0;

    forAll(patchIDs_, i)
    {
        const polyPatch& pp = patches[patchIDs_[i]];

        forAll(pp, facei)
        {
            const face& f = pp.localFaces()[facei];

            vertLabels.append(f.size());
            vtkWriteOps::insert(f + offset, vertLabels);
        }
        offset += pp.nPoints();
    }
    vtkWriteOps::write(os_, binary_, vertLabels);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::patchWriter::writePatchIDs()
{
    const fvMesh& mesh = vMesh_.mesh();

    DynamicList<floatScalar> fField(nFaces_);

    os_ << "patchID 1 " << nFaces_ << " float" << std::endl;

    forAll(patchIDs_, i)
    {
        label patchi = patchIDs_[i];

        const polyPatch& pp = mesh.boundaryMesh()[patchi];

        if (!isA<emptyPolyPatch>(pp))
        {
            vtkWriteOps::insert(scalarField(pp.size(), patchi), fField);
        }
    }
    vtkWriteOps::write(os_, binary_, fField);
}


// ************************************************************************* //
