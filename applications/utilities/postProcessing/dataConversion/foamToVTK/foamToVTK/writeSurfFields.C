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

#include "writeSurfFields.H"
#include "OFstream.H"
#include "floatScalar.H"
#include "vtkWriteFieldOps.H"
#include "emptyFvsPatchFields.H"
#include "fvsPatchFields.H"

// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

void Foam::writeSurfFields
(
    const bool binary,
    const vtkMesh& vMesh,
    const fileName& fileName,
    const UPtrList<const surfaceVectorField>& surfVectorFields
)
{
    const fvMesh& mesh = vMesh.mesh();

    std::ofstream str(fileName.c_str());

    vtkWriteOps::writeHeader
    (
        str,
        binary,
        "surfaceFields"
    );

    str << "DATASET POLYDATA" << std::endl;

    const pointField& fc = mesh.faceCentres();

    str << "POINTS " << mesh.nFaces() << " float" << std::endl;

    DynamicList<floatScalar> pField(3*mesh.nFaces());

    for (label facei = 0; facei < mesh.nFaces(); facei++)
    {
        vtkWriteOps::insert(fc[facei], pField);
    }

    vtkWriteOps::write(str, binary, pField);

    str << "POINT_DATA " << mesh.nFaces() << std::endl
        << "FIELD attributes " << surfVectorFields.size() << std::endl;

    // surfVectorFields
    forAll(surfVectorFields, fieldi)
    {
        const surfaceVectorField& svf = surfVectorFields[fieldi];

        str << svf.name() << " 3 "
            << mesh.nFaces() << " float" << std::endl;

        DynamicList<floatScalar> fField(3*mesh.nFaces());

        for (label facei = 0; facei < mesh.nInternalFaces(); facei++)
        {
            vtkWriteOps::insert(svf[facei], fField);
        }

        forAll(svf.boundaryField(), patchi)
        {
            const fvsPatchVectorField& pf = svf.boundaryField()[patchi];

            const fvPatch& pp = mesh.boundary()[patchi];

            if (isA<emptyFvsPatchVectorField>(pf))
            {
                // Note: loop over polypatch size, not fvpatch size.
                forAll(pp.patch(), i)
                {
                    vtkWriteOps::insert(vector::zero, fField);
                }
            }
            else
            {
                forAll(pf, i)
                {
                    vtkWriteOps::insert(pf[i], fField);
                }
            }
        }

        vtkWriteOps::write(str, binary, fField);
    }
}


// ************************************************************************* //
