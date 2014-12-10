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

#include "processorMeshes.H"
#include "Time.H"
#include "primitiveMesh.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::processorMeshes::read()
{
    // Make sure to clear (and hence unregister) any previously loaded meshes
    // and fields
    forAll(databases_, procI)
    {
        meshes_.set(procI, NULL);
        pointProcAddressing_.set(procI, NULL);
        faceProcAddressing_.set(procI, NULL);
        cellProcAddressing_.set(procI, NULL);
        boundaryProcAddressing_.set(procI, NULL);
    }

    forAll(databases_, procI)
    {
        meshes_.set
        (
            procI,
            new fvMesh
            (
                IOobject
                (
                    meshName_,
                    databases_[procI].timeName(),
                    databases_[procI]
                )
            )
        );

        pointProcAddressing_.set
        (
            procI,
            new labelIOList
            (
                IOobject
                (
                    "pointProcAddressing",
                    meshes_[procI].facesInstance(),
                    meshes_[procI].meshSubDir,
                    meshes_[procI],
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                )
            )
        );

        faceProcAddressing_.set
        (
            procI,
            new labelIOList
            (
                IOobject
                (
                    "faceProcAddressing",
                    meshes_[procI].facesInstance(),
                    meshes_[procI].meshSubDir,
                    meshes_[procI],
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                )
            )
        );

        cellProcAddressing_.set
        (
            procI,
            new labelIOList
            (
                IOobject
                (
                    "cellProcAddressing",
                    meshes_[procI].facesInstance(),
                    meshes_[procI].meshSubDir,
                    meshes_[procI],
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                )
            )
        );

        boundaryProcAddressing_.set
        (
            procI,
            new labelIOList
            (
                IOobject
                (
                    "boundaryProcAddressing",
                    meshes_[procI].facesInstance(),
                    meshes_[procI].meshSubDir,
                    meshes_[procI],
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                )
            )
        );
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::processorMeshes::processorMeshes
(
    PtrList<Time>& databases,
    const word& meshName
)
:
    meshName_(meshName),
    databases_(databases),
    meshes_(databases.size()),
    pointProcAddressing_(databases.size()),
    faceProcAddressing_(databases.size()),
    cellProcAddressing_(databases.size()),
    boundaryProcAddressing_(databases.size())
{
    read();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::fvMesh::readUpdateState Foam::processorMeshes::readUpdate()
{
    fvMesh::readUpdateState stat = fvMesh::UNCHANGED;

    forAll(databases_, procI)
    {
        // Check if any new meshes need to be read.
        fvMesh::readUpdateState procStat = meshes_[procI].readUpdate();

        /*
        if (procStat != fvMesh::UNCHANGED)
        {
            Info<< "Processor " << procI
                << " at time " << databases_[procI].timeName()
                << " detected mesh change " << procStat
                << endl;
        }
        */

        // Combine into overall mesh change status
        if (stat == fvMesh::UNCHANGED)
        {
            stat = procStat;
        }
        else if (stat != procStat)
        {
            FatalErrorIn("processorMeshes::readUpdate()")
                << "Processor " << procI
                << " has a different polyMesh at time "
                << databases_[procI].timeName()
                << " compared to any previous processors." << nl
                << "Please check time " << databases_[procI].timeName()
                << " directories on all processors for consistent"
                << " mesh files."
                << exit(FatalError);
        }
    }

    if
    (
        stat == fvMesh::TOPO_CHANGE
     || stat == fvMesh::TOPO_PATCH_CHANGE
    )
    {
        // Reread all meshes and addresssing
        read();
    }
    return stat;
}


void Foam::processorMeshes::reconstructPoints(fvMesh& mesh)
{
    // Read the field for all the processors
    PtrList<pointIOField> procsPoints(meshes_.size());

    forAll(meshes_, procI)
    {
        procsPoints.set
        (
            procI,
            new pointIOField
            (
                IOobject
                (
                    "points",
                    meshes_[procI].time().timeName(),
                    polyMesh::meshSubDir,
                    meshes_[procI],
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                )
            )
        );
    }

    // Create the new points
    vectorField newPoints(mesh.nPoints());

    forAll(meshes_, procI)
    {
        const vectorField& procPoints = procsPoints[procI];

        // Set the cell values in the reconstructed field

        const labelList& pointProcAddressingI = pointProcAddressing_[procI];

        if (pointProcAddressingI.size() != procPoints.size())
        {
            FatalErrorIn("processorMeshes")
                << "problem :"
                << " pointProcAddressingI:" << pointProcAddressingI.size()
                << " procPoints:" << procPoints.size()
                << abort(FatalError);
        }

        forAll(pointProcAddressingI, pointI)
        {
            newPoints[pointProcAddressingI[pointI]] = procPoints[pointI];
        }
    }

    mesh.movePoints(newPoints);
    mesh.write();
}


// ************************************************************************* //
