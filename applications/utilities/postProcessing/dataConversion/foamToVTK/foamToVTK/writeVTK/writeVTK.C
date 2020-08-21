/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2020 OpenFOAM Foundation
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

#include "writeVTK.H"
#include "Time.H"
#include "vtkMesh.H"
#include "internalWriter.H"
#include "OSspecific.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(writeVTK, 0);
    addToRunTimeSelectionTable(functionObject, writeVTK, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::writeVTK::writeVTK
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    objectNames_()
{
    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::writeVTK::~writeVTK()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::writeVTK::read(const dictionary& dict)
{
    dict.lookup("objects") >> objectNames_;

    return true;
}


bool Foam::functionObjects::writeVTK::execute()
{
    return true;
}


bool Foam::functionObjects::writeVTK::write()
{
    Info<< type() << " " << name() << " output:" << nl;

    Info<< "Time: " << time_.timeName() << endl;

    word timeDesc = time_.timeName();

    // VTK/ directory in the case
    fileName fvPath(time_.path()/"VTK");

    mkDir(fvPath);

    string vtkName = time_.caseName();

    if (Pstream::parRun())
    {
        // Strip off leading casename, leaving just processor_DDD ending.
        string::size_type i = vtkName.rfind("processor");

        if (i != string::npos)
        {
            vtkName = vtkName.substr(i);
        }
    }

    // Create file and write header
    fileName vtkFileName
    (
        fvPath/vtkName
      + "_"
      + timeDesc
      + ".vtk"
    );

    Info<< "    Internal  : " << vtkFileName << endl;

    vtkMesh vMesh(const_cast<fvMesh&>(mesh_));

    // Write mesh
    internalWriter writer(vMesh, false, vtkFileName);

    UPtrList<const volScalarField> vsf(lookupFields<volScalarField>());
    UPtrList<const volVectorField> vvf(lookupFields<volVectorField>());
    UPtrList<const volSphericalTensorField> vsptf
    (
        lookupFields<volSphericalTensorField>()
    );
    UPtrList<const volSymmTensorField> vstf(lookupFields<volSymmTensorField>());
    UPtrList<const volTensorField> vtf(lookupFields<volTensorField>());

    // Write header for cellID and volFields
    vtkWriteOps::writeCellDataHeader
    (
        writer.os(),
        vMesh.nFieldCells(),
        1 + vsf.size() + vvf.size() + vsptf.size() + vstf.size() + vtf.size()
    );

    // Write cellID field
    writer.writeCellIDs();

    // Write volFields
    writer.write(vsf);
    writer.write(vvf);
    writer.write(vsptf);
    writer.write(vstf);
    writer.write(vtf);

    return true;
}


// ************************************************************************* //
