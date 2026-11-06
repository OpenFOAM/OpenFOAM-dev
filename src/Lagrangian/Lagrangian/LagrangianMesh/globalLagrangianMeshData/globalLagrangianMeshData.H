/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2026 OpenFOAM Foundation
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

Class
    Foam::globalLagrangianMeshData

Description
    Global Lagrangian mesh data. Minimal wrapper just to provide interface
    consistency for the fields.

SourceFiles
    globalLagrangianMeshData.C

\*---------------------------------------------------------------------------*/

#ifndef globalLagrangianMeshData_H
#define globalLagrangianMeshData_H

#include "className.H"
#include "lduSchedule.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class LagrangianMesh;

/*---------------------------------------------------------------------------*\
                     Class globalLagrangianMeshData Declaration
\*---------------------------------------------------------------------------*/

class globalLagrangianMeshData
{
    // Private Data

        //- Reference to the mesh
        const LagrangianMesh& mesh_;


public:

    //- Runtime type information
    ClassName("globalLagrangianMeshData");


    // Constructors

        //- Construct for a mesh
        globalLagrangianMeshData(const LagrangianMesh& mesh);

        //- Disallow default bitwise copy construction
        globalLagrangianMeshData(const globalLagrangianMeshData&) = delete;


    //- Destructor
    ~globalLagrangianMeshData();


    // Member Functions

        //- Use the same patch schedule as the polyhedral mesh
        const lduSchedule& patchSchedule() const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const globalLagrangianMeshData&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
