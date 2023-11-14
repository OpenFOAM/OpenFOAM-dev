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
    Foam::pointFieldDecomposer

Description
    Point field decomposer.

SourceFiles
    pointFieldDecomposer.C
    pointFieldDecomposerTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef pointFieldDecomposer_H
#define pointFieldDecomposer_H

#include "pointMesh.H"
#include "pointFields.H"
#include "forwardFieldMapper.H"
#include "IOobjectList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class fvMesh;

/*---------------------------------------------------------------------------*\
                    Class pointFieldDecomposer Declaration
\*---------------------------------------------------------------------------*/

class pointFieldDecomposer
{
public:

    // Public classes

        //- Point patch field decomposer class
        class patchFieldDecomposer
        :
            public labelList,
            public forwardFieldMapper
        {
            // Private Member Functions

                //- Generate the addressing
                static labelList addressing
                (
                    const pointPatch& completePatch,
                    const pointPatch& procPatch,
                    const labelList& pointProcAddressing
                );


        public:

            // Constructors

                //- Construct given patches and addressing
                patchFieldDecomposer
                (
                    const pointPatch& completePatch,
                    const pointPatch& procPatch,
                    const labelList& pointProcAddressing
                );
        };


private:

    // Private Data

        //- Reference to complete mesh
        const pointMesh& completeMesh_;

        //- Reference to processor meshes
        const PtrList<fvMesh>& procMeshes_;

        //- Reference to point addressing
        const labelListList& pointProcAddressing_;

        //- List of patch field decomposers
        PtrList<PtrList<patchFieldDecomposer>> patchFieldDecomposers_;


    // Private Member Functions

        //- Decompose a field
        template<class Type>
        PtrList<PointField<Type>>
        decomposeField(const IOobject& fieldIoObject) const;


public:

    // Constructors

        //- Construct from components
        pointFieldDecomposer
        (
            const pointMesh& completeMesh,
            const PtrList<fvMesh>& procMeshes,
            const labelListList& pointAddressing
        );

        //- Disallow default bitwise copy construction
        pointFieldDecomposer(const pointFieldDecomposer&) = delete;


    //- Destructor
    ~pointFieldDecomposer();


    // Member Functions

        //- Return whether anything in the object list gets decomposed
        static bool decomposes(const IOobjectList& objects);

        //- Read, decompose and write all fields
        template<class Type>
        void decomposeFields(const IOobjectList& objects);


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const pointFieldDecomposer&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "pointFieldDecomposerTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
