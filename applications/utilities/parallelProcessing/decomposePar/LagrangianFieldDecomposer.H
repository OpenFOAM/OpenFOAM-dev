/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2025 OpenFOAM Foundation
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
    Foam::LagrangianFieldDecomposer

Description
    Lagrangian field decomposer

SourceFiles
    LagrangianFieldDecomposer.C
    LagrangianFieldDecomposerTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef LagrangianFieldDecomposer_H
#define LagrangianFieldDecomposer_H

#include "fvMesh.H"
#include "LagrangianMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class IOobjectList;

/*---------------------------------------------------------------------------*\
                  Class LagrangianFieldDecomposer Declaration
\*---------------------------------------------------------------------------*/

class LagrangianFieldDecomposer
{
    // Private Data

        //- The complete mesh
        LagrangianMesh completeMesh_;

        //- List of processor meshes
        PtrList<LagrangianMesh> procMeshes_;

        //- For each processor, for each element, the global element index
        labelListList particleProcAddressing_;


    // Private Member Functions

        //- Return whether anything in the object list gets decomposed
        template<class GeoField>
        static bool decomposes(const IOobjectList& objects);

        //- Return the decomposed internal field given a list of processor
        //  internal fields
        template<class Type, template<class> class PrimitiveField>
        PtrList<DimensionedField<Type, LagrangianMesh, PrimitiveField>>
        decomposeLagrangianField
        (
            const DimensionedField<Type, LagrangianMesh, PrimitiveField>&
                completeField
        ) const;

        //- Return the decomposed field given a list of processor fields
        template<class Type, template<class> class PrimitiveField>
        PtrList<GeometricField<Type, LagrangianMesh, PrimitiveField>>
        decomposeLagrangianField
        (
            const GeometricField<Type, LagrangianMesh, PrimitiveField>&
                completeField
        ) const;

        //- decompose a field
        template<class GeoField>
        PtrList<GeoField> decomposeLagrangianField(const IOobject&) const;


public:

    // Constructors

        //- Construct from components
        LagrangianFieldDecomposer
        (
            const fvMesh& completeFvMesh,
            const PtrList<fvMesh>& procFvMeshes,
            const labelListList& faceProcAddressing,
            const labelListList& cellProcAddressing,
            const word& LagrangianName
        );

        //- Disallow default bitwise copy construction
        LagrangianFieldDecomposer
        (
            const LagrangianFieldDecomposer&
        ) = delete;


    //- Destructor
    ~LagrangianFieldDecomposer();


    // Member Functions

        //- Return whether anything in the object list gets decomposed
        static bool decomposes(const IOobjectList& objects);

        //- Read, decompose and write all Lagrangian fields
        template<class GeoField>
        void decomposeFields(const IOobjectList& objects) const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const LagrangianFieldDecomposer&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "LagrangianFieldDecomposerTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
