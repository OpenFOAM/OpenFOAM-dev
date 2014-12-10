/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013-2014 OpenFOAM Foundation
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
    createExternalCoupledPatchGeometry.

Description
    Application to generate the patch geometry (points and faces) for use
    with the externalCoupled boundary condition.

    Usage:

        createExternalCoupledPatchGeometry \<fieldName\>

    On execution, the field \<fieldName\> is read, and its boundary conditions
    interrogated for the presence of an \c externalCoupled type.  If found,
    the patch geometry (points and faces) for the coupled patches are output
    to the communications directory.

Note:
    The addressing is patch-local, i.e. point indices for each patch point
    used for face addressing starts at index 0.

SeeAlso
    externalCoupledMixedFvPatchField

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "createExternalCoupledPatchGeometryTemplates.H"
#include "IOobjectList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "addRegionOption.H"
    argList::validArgs.append("fieldName");
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createNamedMesh.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    const word fieldName = args[1];

    IOobjectList objects(IOobjectList(mesh, mesh.time().timeName()));

    label processed = -1;
    processField<scalar>(mesh, objects, fieldName, processed);
    processField<vector>(mesh, objects, fieldName, processed);
    processField<sphericalTensor>(mesh, objects, fieldName, processed);
    processField<symmTensor>(mesh, objects, fieldName, processed);
    processField<tensor>(mesh, objects, fieldName, processed);

    if (processed == -1)
    {
        Info<< "Field " << fieldName << " not found" << endl;
    }

    Info<< "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
