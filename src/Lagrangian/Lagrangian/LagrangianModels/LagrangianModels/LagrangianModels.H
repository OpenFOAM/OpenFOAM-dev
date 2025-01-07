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
    Foam::LagrangianModels

Description
    List of Lagrangian models, constructed as a (Lagrangian) mesh object.
    Provides similar functions to the models themselves and forwards them to
    each model in turn. This is the high level model interface used by clouds
    when constructing their injections and transport equations.

SourceFiles
    LagrangianModels.C
    LagrangianModelsTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef LagrangianModels_H
#define LagrangianModels_H

#include "LagrangianModel.H"
#include "PtrListDictionary.H"
#include "DemandDrivenMeshObject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class LagrangianModels Declaration
\*---------------------------------------------------------------------------*/

class LagrangianModels
:
    public DemandDrivenMeshObject
    <
        LagrangianMesh,
        TopoChangeableMeshObject,
        LagrangianModels
    >,
    public dictionary,
    public PtrListDictionary<LagrangianModel>
{
    // Private Member Data

        //- Time index to check that all defined sources have been applied
        mutable label checkTimeIndex_;

        //- Sets of the fields that have had sources added by the
        //  LagrangianModels
        mutable PtrList<wordHashSet> addSupFields_;


    // Private Member Functions

        //- Create IO object for an optional LagrangianModels file
        IOobject io(const LagrangianMesh& mesh) const;

        //- Helper for modelTypeFieldSourceTypes. Inserts the field source type
        //  into to the table if the model is of the required type.
        template<class ... ModelAndFieldSourceTypes>
        struct modelTypeFieldSourceType;

        //- Check that all LagrangianModels have been applied
        void checkApplied() const;

        //- Return a source for an equation
        template
        <
            class Type,
            template<class> class PrimitiveField,
            class ... AlphaRhoFieldTypes
        >
        tmp<LagrangianEqn<Type>> sourceTerm
        (
            const LagrangianSubField<Type, PrimitiveField>& eqnField,
            const LagrangianSubScalarField& deltaT,
            const AlphaRhoFieldTypes& ... alphaRhoFields
        ) const;


protected:

    //- Let the mesh object use the protected constructors
    friend class DemandDrivenMeshObject
    <
        LagrangianMesh,
        TopoChangeableMeshObject,
        LagrangianModels
    >;


    // Protected Constructors

        //- Construct for a mesh
        explicit LagrangianModels(const LagrangianMesh& mesh);


public:

    //- Runtime type information
    TypeName("LagrangianModels");


    // Constructors

        //- Disallow default bitwise copy construction
        LagrangianModels(const LagrangianModels&) = delete;

        //- Inherit the base New method
        using DemandDrivenMeshObject
        <
            LagrangianMesh,
            TopoChangeableMeshObject,
            LagrangianModels
        >::New;


    //- Destructor
    virtual ~LagrangianModels();


    // Member Functions

        //- Declare LagrangianModels to be a global dictionary
        virtual bool global() const
        {
            return true;
        }


        // Checks

            //- Return true if the LagrangianModels adds a source term to the
            //  given field's transport equation
            bool addsSupToField(const word&) const;

            //- Return true if the LagrangianModels adds a source term to the
            //  given field's transport equation
            template<class Type, template<class> class PrimitiveField>
            bool addsSupToField
            (
                const LagrangianSubField<Type, PrimitiveField>& field
            ) const;


        //- Correct the LagrangianModels
        void correct();

        //- Identify elements in the Lagrangian mesh which are to be
        //  instantaneously modified and put them in contiguous groups
        LagrangianSubMesh preModify(LagrangianMesh& mesh) const;

        //- Instantaneously modify and/or create and remove elements in the
        //  Lagrangian mesh
        LagrangianSubMesh modify
        (
            LagrangianMesh& mesh,
            const LagrangianSubMesh& modifiedMesh
        ) const;

        //- Solve equations and/or update continually changing properties
        void calculate
        (
            const LagrangianSubScalarField& deltaT,
            const bool final
        );

        //- Return a table of field source types that are chosen to match given
        //  model types. So, e.g., a zero field source for every injection
        //  model. Useful for programmatically constructing fields. Template
        //  arguments should alternate between model and field types; i.e.,
        //  model-type-1, field-source-type-1, model-type-2,
        //  field-source-type-2, model-type-3, ...
        template<class ... ModelAndFieldSourceTypes>
        HashTable<word> modelTypeFieldSourceTypes() const;


        // Sources

            //- Return the fractional source
            tmp<LagrangianSubScalarField> source
            (
                const LagrangianSubScalarField& deltaT
            ) const;

            //- Return source for an equation
            template<class Type, template<class> class PrimitiveField>
            tmp<LagrangianEqn<Type>> source
            (
                const LagrangianSubScalarField& deltaT,
                const LagrangianSubField<Type, PrimitiveField>& field
            ) const;

            //- Return source for an equation
            template
            <
                class Type,
                template<class> class PrimitiveField,
                template<class> class PrimitiveEqnField
            >
            tmp<LagrangianEqn<Type>> sourceProxy
            (
                const LagrangianSubScalarField& deltaT,
                const LagrangianSubField<Type, PrimitiveField>& field,
                const LagrangianSubField<Type, PrimitiveEqnField>& eqnField
            ) const;

            //- Return source for a mass-weighted equation
            template<class Type, template<class> class PrimitiveField>
            tmp<LagrangianEqn<Type>> source
            (
                const LagrangianSubScalarField& deltaT,
                const LagrangianSubField<scalar, PrimitiveField>& m,
                const LagrangianSubField<Type, PrimitiveField>& field
            ) const;

            //- Return source for a mass-weighted equation
            template
            <
                class Type,
                template<class> class PrimitiveField,
                template<class> class PrimitiveEqnField
            >
            tmp<LagrangianEqn<Type>> sourceProxy
            (
                const LagrangianSubScalarField& deltaT,
                const LagrangianSubField<scalar, PrimitiveField>& m,
                const LagrangianSubField<Type, PrimitiveField>& field,
                const LagrangianSubField<Type, PrimitiveEqnField>& eqnField
            ) const;


        // Mesh changes

            //- Update for mesh motion. Only for mesh object. Does nothing.
            //  Lagrangian evolves continuously across the range of mesh
            //  motion, so no instantaneous update is needed.
            virtual bool movePoints();

            //- Update topology using the given map
            virtual void topoChange(const polyTopoChangeMap&);

            //- Update from another mesh using the given map
            virtual void mapMesh(const polyMeshMap&);

            //- Redistribute or update using the given distribution map
            virtual void distribute(const polyDistributionMap&);


        // IO

            //- Read the LagrangianModels dictionary if it has changed
            //  and update the models
            virtual bool read();

            //- ReadData function which reads the LagrangianModels dictionary
            //  required for regIOobject read operation
            virtual bool readData(Istream&);

            //- writeData function required by regIOobject but not used
            //  for this class, writeObject is used instead
            virtual bool writeData(Ostream& os) const;

            //- Write the LagrangianModels
            virtual bool writeObject
            (
                IOstream::streamFormat fmt,
                IOstream::versionNumber ver,
                IOstream::compressionType cmp,
                const bool write
            ) const;


    // Member Operators

        //- Inherit the PtrListDictionary index operators
        using PtrListDictionary<LagrangianModel>::operator[];

        //- Inherit the PtrListDictionary size function
        using PtrListDictionary<LagrangianModel>::size;

        //- Disallow default bitwise assignment
        void operator=(const LagrangianModels&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- Trait for obtaining global status
template<>
struct typeGlobal<LagrangianModels>
{
    static const bool global = true;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "LagrangianModelsTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
