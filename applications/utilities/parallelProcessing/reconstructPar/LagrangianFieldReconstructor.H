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
    Foam::LagrangianFieldReconstructor

Description
    Lagrangian field reconstructor

SourceFiles
    LagrangianFieldReconstructor.C
    LagrangianFieldReconstructorTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef LagrangianFieldReconstructor_H
#define LagrangianFieldReconstructor_H

#include "fvMesh.H"
#include "LagrangianMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class IOobjectList;

/*---------------------------------------------------------------------------*\
                  Class LagrangianFieldReconstructor Declaration
\*---------------------------------------------------------------------------*/

class LagrangianFieldReconstructor
{
    // Private Data

        //- The complete mesh
        LagrangianMesh completeMesh_;

        //- List of processor meshes
        PtrList<LagrangianMesh> procMeshes_;


    // Private Member Functions

        //- Return whether anything in the object list gets reconstructed
        template<class GeoField>
        static bool reconstructs
        (
            const IOobjectList& objects,
            const HashSet<word>& selectedFields
        );

        //- Return the reconstructed primitive field given a list of processor
        //  (internal) fields
        template
        <
            class Type,
            template<class> class PrimitiveField,
            template<class, class, template<class> class> class GeoField
        >
        tmp<PrimitiveField<Type>> reconstructLagrangianPrimitiveField
        (
            const PtrList
            <
                GeoField<Type, LagrangianMesh, PrimitiveField>
            >& procFields
        ) const;

        //- Return the reconstructed internal field given a list of processor
        //  internal fields
        template<class Type, template<class> class PrimitiveField>
        tmp<DimensionedField<Type, LagrangianMesh, PrimitiveField>>
        reconstructLagrangianField
        (
            const PtrList
            <
                DimensionedField<Type, LagrangianMesh, PrimitiveField>
            >& procFields
        ) const;

        //- Return the reconstructed field given a list of processor fields
        template<class Type, template<class> class PrimitiveField>
        tmp<GeometricField<Type, LagrangianMesh, PrimitiveField>>
        reconstructLagrangianField
        (
            const PtrList
            <
                GeometricField<Type, LagrangianMesh, PrimitiveField>
            >& procFields
        ) const;

        //- Reconstruct a field
        template<class GeoField>
        tmp<GeoField> reconstructLagrangianField(const IOobject&) const;


public:

    // Constructors

        //- Construct from components
        LagrangianFieldReconstructor
        (
            const fvMesh& completeFvMesh,
            const PtrList<fvMesh>& procFvMeshes,
            const labelListList& faceProcAddressing,
            const labelListList& cellProcAddressing,
            const word& LagrangianName
        );

        //- Disallow default bitwise copy construction
        LagrangianFieldReconstructor
        (
            const LagrangianFieldReconstructor&
        ) = delete;


    //- Destructor
    ~LagrangianFieldReconstructor();


    // Member Functions

        //- Return whether anything in the object list gets reconstructed
        static bool reconstructs
        (
            const IOobjectList& objects,
            const HashSet<word>& selectedFields
        );

        //- Read, reconstruct and write all/selected Lagrangian fields
        template<class GeoField>
        void reconstructFields
        (
            const IOobjectList& objects,
            const HashSet<word>& selectedFields
        ) const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const LagrangianFieldReconstructor&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "LagrangianFieldReconstructorTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
