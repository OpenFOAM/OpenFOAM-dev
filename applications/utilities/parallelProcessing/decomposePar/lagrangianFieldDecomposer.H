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
    Foam::lagrangianFieldDecomposer

Description
    Lagrangian field decomposer.

SourceFiles
    lagrangianFieldDecomposer.C
    lagrangianFieldDecomposerTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef lagrangianFieldDecomposer_H
#define lagrangianFieldDecomposer_H

#include "cloud.H"
#include "fvMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class IOobjectList;

/*---------------------------------------------------------------------------*\
                  Class lagrangianFieldDecomposer Declaration
\*---------------------------------------------------------------------------*/

class lagrangianFieldDecomposer
{
    // Private Data

        //- Reference to complete mesh
        const fvMesh& completeMesh_;

        //- List of processor meshes
        const PtrList<fvMesh>& procMeshes_;

        //- For each processor, for each particle, the global particle index
        labelListList particleProcAddressing_;

        //- The name of the cloud
        const word cloudName_;


    // Private Member Functions

        //- Decompose a field
        template<class Type, template<class> class IOContainer>
        PtrList<IOContainer<Type>>
        decomposeField(const IOobject& fieldIoObject) const;

        //- Read, decompose and write all fields
        template
        <
            class Type,
            template<class> class IOContainer,
            template<class> class IOContainerType
        >
        void decomposeFields(const IOobjectList& objects);


public:

    // Constructors

        //- Construct from components
        lagrangianFieldDecomposer
        (
            const fvMesh& completeMesh,
            const PtrList<fvMesh>& procMeshes,
            const labelListList& faceProcAddressing,
            const labelListList& cellProcAddressing,
            const word& cloudName
        );

        //- Disallow default bitwise copy construction
        lagrangianFieldDecomposer(const lagrangianFieldDecomposer&) = delete;


    //- Destructor
    ~lagrangianFieldDecomposer();


    // Member Functions

        //- Return whether anything in the object list gets decomposed
        static bool decomposes(const IOobjectList& objects);

        //- Read, decompose and write all fields
        template<class Type>
        void decomposeFields(const IOobjectList& objects);


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const lagrangianFieldDecomposer&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "lagrangianFieldDecomposerTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
