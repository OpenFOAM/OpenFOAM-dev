/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2023 OpenFOAM Foundation
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
    Foam::lagrangianFieldReconstructor

Description
    Lagrangian field reconstructor.

SourceFiles
    lagrangianFieldReconstructor.C
    lagrangianFieldReconstructorTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef lagrangianFieldReconstructor_H
#define lagrangianFieldReconstructor_H

#include "cloud.H"
#include "fvMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class IOobjectList;

/*---------------------------------------------------------------------------*\
                  Class lagrangianFieldReconstructor Declaration
\*---------------------------------------------------------------------------*/

class lagrangianFieldReconstructor
{
    // Private Data

        //- Reference to complete mesh
        const fvMesh& completeMesh_;

        //- List of processor meshes
        const PtrList<fvMesh>& procMeshes_;

        //- The name of the cloud
        const word cloudName_;


    // Private Member Functions

        //- Return whether anything in the object list gets reconstructed
        template<class FieldType>
        static bool reconstructs
        (
            const IOobjectList& objects,
            const HashSet<word>& selectedFields
        );

        //- Reconstruct a field
        template
        <
            class Type,
            template<class> class IOContainer,
            template<class> class IOContainerType
        >
        tmp<IOContainer<Type>>
        reconstructField(const IOobject& fieldIoObject) const;

        //- Read, reconstruct and write all fields
        template
        <
            class Type,
            template<class> class IOContainer,
            template<class> class IOContainerType
        >
        void reconstructFields
        (
            const IOobjectList& objects,
            const HashSet<word>& selectedFields
        );


public:

    // Constructors

        //- Construct from components
        lagrangianFieldReconstructor
        (
            const fvMesh& completeMesh,
            const PtrList<fvMesh>& procMeshes,
            const labelListList& faceProcAddressing,
            const labelListList& cellProcAddressing,
            const word& cloudName
        );

        //- Disallow default bitwise copy construction
        lagrangianFieldReconstructor
        (
            const lagrangianFieldReconstructor&
        ) = delete;


    //- Destructor
    ~lagrangianFieldReconstructor();


    // Member Functions

        //- Return whether anything in the object list gets reconstructed
        static bool reconstructs
        (
            const IOobjectList& objects,
            const HashSet<word>& selectedFields
        );

        //- Read, reconstruct and write all fields
        template<class Type>
        void reconstructFields
        (
            const IOobjectList& objects,
            const HashSet<word>& selectedFields
        );


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const lagrangianFieldReconstructor&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "lagrangianFieldReconstructorTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
