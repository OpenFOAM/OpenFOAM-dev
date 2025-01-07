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
    Foam::LagrangianSubMesh

Description
    Mesh that relates to a sub-section of a Lagrangian mesh. This is used to
    construct fields that relate to a contiguous sub-set of the Lagrangian
    elements. This class only stores references and the indices defining the
    range of the sub-set, so it is very lightweight and can be constructed and
    thrown away largely without consideration of expense.

SourceFiles
    LagrangianSubMesh.C

\*---------------------------------------------------------------------------*/

#ifndef LagrangianSubMesh_H
#define LagrangianSubMesh_H

#include "GeoMesh.H"
#include "LagrangianState.H"
#include "LagrangianFieldsFwd.H"
#include "LagrangianSubFieldsFwd.H"
#include "polyMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class LagrangianMesh;

/*---------------------------------------------------------------------------*\
                         Class LagrangianSubMesh Declaration
\*---------------------------------------------------------------------------*/

class LagrangianSubMesh
:
    public GeoMesh<polyMesh>
{
    // Private Data

        //- Reference to the Lagrangian mesh
        const LagrangianMesh& mesh_;

        //- Group of the elements in this subset
        const LagrangianGroup group_;

        //- Size of the subset
        label size_;

        //- Start index of the subset
        label start_;

        //- Index with which to compare instances
        const label index_;


    // Private Constructors

        //- Construct from components
        explicit LagrangianSubMesh
        (
            const LagrangianMesh& mesh,
            const LagrangianGroup group,
            const label size,
            const label start,
            const label index
        );


public:

    // Friend classes

        //- Allow the Lagrangian mesh to construct with a specified index
        friend class LagrangianMesh;


    //- Runtime type information
    ClassName("LagrangianSubMesh");


    // Public Type Definitions

        //- Mesh type
        typedef LagrangianSubMesh Mesh;


    // Constructors

        //- Construct from components, except for the index which is
        //  automatically generated from the complete mesh
        explicit LagrangianSubMesh
        (
            const LagrangianMesh& mesh,
            const LagrangianGroup group,
            const label size,
            const label start
        );

        //- Construct from Lagrangian mesh, group offsets and group
        explicit LagrangianSubMesh
        (
            const LagrangianMesh& mesh,
            const labelList& groupOffsets,
            const LagrangianGroup group
        );


    //- Destructor
    ~LagrangianSubMesh();


    // Member functions

        // Access

            //- Return the mesh
            inline const LagrangianMesh& mesh() const;

            //- Return the group
            inline LagrangianGroup group() const;

            //- Return size
            inline label size() const;

            //- Return size
            inline label globalSize() const;

            //- Return size
            static inline label size(const LagrangianSubMesh& subMesh);

            //- Return whether or not the mesh is empty
            inline bool empty() const;

            //- Return start
            inline label start() const;

            //- Return end
            inline label end() const;

            //- Return the index
            inline label index() const;


        // Sub-setting

            //- Return a sub-list corresponding to this sub-mesh
            template<class Type>
            SubList<Type> sub(const List<Type>& list) const;

            //- Return a sub-field corresponding to this sub-mesh
            template<class Type>
            SubField<Type> sub(const Field<Type>& field) const;

            //- Return a sub-dimensioned-field corresponding to this sub-mesh
            template<class Type, template<class> class PrimitiveField>
            tmp<LagrangianSubSubField<Type>> sub
            (
                const DimensionedField<Type, LagrangianMesh, PrimitiveField>&
            ) const;


        // Geometry

            //- Return the face normals at the Lagrangian locations
            template<class FieldType>
            static tmp<FieldType> nf
            (
                const LagrangianSubScalarSubField& fraction
            );

            //- Return the face normals at the Lagrangian locations
            template<class FieldType>
            tmp<FieldType> nf
            (
                const LagrangianScalarInternalDynamicField& fraction
            ) const;

            //- Return the face velocities at the Lagrangian locations
            template<class FieldType>
            static tmp<FieldType> Uf
            (
                const LagrangianSubScalarSubField& fraction
            );

            //- Return the face velocities at the Lagrangian locations
            template<class FieldType>
            tmp<FieldType> Uf
            (
                const LagrangianScalarInternalDynamicField& fraction
            ) const;


    // Member Operators

        //- Add a sub-mesh to this one. Must relate to adjacent elements.
        void operator+=(const LagrangianSubMesh&);
};


// Template Specialisations

template<>
tmp<vectorField> LagrangianSubMesh::nf
(
    const LagrangianSubScalarSubField& fraction
);

template<>
tmp<LagrangianSubVectorField> LagrangianSubMesh::nf
(
    const LagrangianSubScalarSubField& fraction
);

template<>
tmp<vectorField> LagrangianSubMesh::nf
(
    const LagrangianScalarInternalDynamicField& fraction
) const;

template<>
tmp<LagrangianSubVectorField> LagrangianSubMesh::nf
(
    const LagrangianScalarInternalDynamicField& fraction
) const;

template<>
tmp<vectorField> LagrangianSubMesh::Uf
(
    const LagrangianSubScalarSubField& fraction
);

template<>
tmp<LagrangianSubVectorField> LagrangianSubMesh::Uf
(
    const LagrangianSubScalarSubField& fraction
);

template<>
tmp<vectorField> LagrangianSubMesh::Uf
(
    const LagrangianScalarInternalDynamicField& fraction
) const;

template<>
tmp<LagrangianSubVectorField> LagrangianSubMesh::Uf
(
    const LagrangianScalarInternalDynamicField& fraction
) const;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "LagrangianSubMeshI.H"

#ifdef NoRepository
    #include "LagrangianSubMeshTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
