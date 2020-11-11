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

#include "lagrangianWriter.H"
#include "vtkWriteFieldOps.H"
#include "Cloud.H"
#include "passiveParticle.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::lagrangianWriter::lagrangianWriter
(
    const vtkMesh& vMesh,
    const bool binary,
    const fileName& fName,
    const word& cloudName,
    const bool dummyCloud
)
:
    vMesh_(vMesh),
    binary_(binary),
    fName_(fName),
    cloudName_(cloudName),
    os_(fName.c_str())
{
    const fvMesh& mesh = vMesh_.mesh();

    // Write header
    vtkWriteOps::writeHeader(os_, binary_, mesh.time().caseName());
    os_ << "DATASET POLYDATA" << std::endl;

    if (dummyCloud)
    {
        nParcels_ = 0;

        os_ << "POINTS " << nParcels_ << " float" << std::endl;

        os_ << "VERTICES " << nParcels_ << ' ' << 2*nParcels_ << std::endl;
    }
    else
    {
        Cloud<passiveParticle> parcels(mesh, cloudName_, false);

        nParcels_ = parcels.size();

        os_ << "POINTS " << nParcels_ << " float" << std::endl;

        DynamicList<floatScalar> partField(3*parcels.size());
        forAllConstIter(Cloud<passiveParticle>, parcels, elmnt)
        {
            vtkWriteOps::insert(elmnt().position(), partField);
        }
        vtkWriteOps::write(os_, binary_, partField);

        os_ << "VERTICES " << nParcels_ << ' ' << 2*nParcels_ << std::endl;

        DynamicList<label> vertexPoints(2*parcels.size());
        forAll(parcels, parceli)
        {
            vertexPoints.append(1);
            vertexPoints.append(parceli);
        }
        vtkWriteOps::write(os_, binary, vertexPoints);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::lagrangianWriter::writeFieldsHeader(const label nFields)
{
    os_ << "POINT_DATA " << nParcels_ << std::endl
        << "FIELD attributes " << nFields
        << std::endl;
}


// ************************************************************************* //
