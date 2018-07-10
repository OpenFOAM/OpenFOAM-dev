/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
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

\*---------------------------------------------------------------------------*/

#include "fvFieldReconstructor.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fvFieldReconstructor::fvFieldReconstructor
(
    fvMesh& mesh,
    const PtrList<fvMesh>& procMeshes,
    const PtrList<labelIOList>& faceProcAddressing,
    const PtrList<labelIOList>& cellProcAddressing,
    const PtrList<labelIOList>& boundaryProcAddressing
)
:
    mesh_(mesh),
    procMeshes_(procMeshes),
    faceProcAddressing_(faceProcAddressing),
    cellProcAddressing_(cellProcAddressing),
    boundaryProcAddressing_(boundaryProcAddressing),
    nReconstructed_(0)
{
    forAll(procMeshes_, proci)
    {
        const fvMesh& procMesh = procMeshes_[proci];
        if
        (
            faceProcAddressing[proci].size() != procMesh.nFaces()
         || cellProcAddressing[proci].size() != procMesh.nCells()
         || boundaryProcAddressing[proci].size() != procMesh.boundary().size()
        )
        {
            FatalErrorInFunction
                << "Size of maps does not correspond to size of mesh"
                << " for processor " << proci << endl
                << "faceProcAddressing : " << faceProcAddressing[proci].size()
                << " nFaces : " << procMesh.nFaces() << endl
                << "cellProcAddressing : " << cellProcAddressing[proci].size()
                << " nCell : " << procMesh.nCells() << endl
                << "boundaryProcAddressing : "
                << boundaryProcAddressing[proci].size()
                << " nFaces : " << procMesh.boundary().size()
                << exit(FatalError);
        }
    }
}


// ************************************************************************* //
