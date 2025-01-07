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
    Foam::LagrangianModel

Description
    Base class for Lagrangian models

SourceFiles
    LagrangianModelI.H
    LagrangianModel.C

\*---------------------------------------------------------------------------*/

#ifndef LagrangianModel_H
#define LagrangianModel_H

#include "LagrangianEqn.H"
#include "LagrangianModelM.H"
#include "stateModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class LagrangianMesh;
class polyTopoChangeMap;
class polyMeshMap;
class polyDistributionMap;

/*---------------------------------------------------------------------------*\
                          Class LagrangianModel Declaration
\*---------------------------------------------------------------------------*/

class LagrangianModel
:
    public stateModel
{
public:

    // Public Enumerations

        //- Enumeration of the types of instantaneous modification
        enum class modification : label
        {
            change = 0,
            remove = 1
        };


    // Public Type Definitions

        //- Class containing an element-index and a modification-enumeration
        typedef LagrangianMesh::elementGroup<modification> elementModification;


private:

    // Private Member Data

        //- Model name
        const word name_;

        //- Reference to the mesh
        const LagrangianMesh& mesh_;


protected:

    // Protected Member Functions

        //- Add a source term to an equation
        template<class Type>
        void addSupType
        (
            const LagrangianSubScalarField& deltaT,
            const LagrangianSubSubField<Type>& field,
            LagrangianEqn<Type>& eqn
        ) const;

        //- Add a source term to a mass-weighted equation
        template<class Type>
        void addSupType
        (
            const LagrangianSubScalarField& deltaT,
            const LagrangianSubScalarSubField& m,
            const LagrangianSubSubField<Type>& field,
            LagrangianEqn<Type>& eqn
        ) const;


public:

    //- Runtime type information
    TypeName("LagrangianModel");


    //- Declare run-time constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        LagrangianModel,
        dictionary,
        (
            const word& name,
            const LagrangianMesh& mesh,
            const dictionary& modelDict,
            const dictionary& stateDict
        ),
        (name, mesh, modelDict, stateDict)
    );


    // Static Member Functions

        //- Return the name of the field associated with a source term
        template<class AlphaRhoFieldType, class ... AlphaRhoFieldTypes>
        static word fieldName
        (
            const AlphaRhoFieldType& alphaRhoField,
            const AlphaRhoFieldTypes& ... alphaRhoFields
        );

        //- Return the name of the field associated with a source term (base
        //  condition for the above)
        template<class AlphaRhoFieldType>
        static word fieldName(const AlphaRhoFieldType& alphaRhoField);

        //- Return the name of the product of the fields associated with a
        //  source term
        template<class AlphaRhoFieldType, class ... AlphaRhoFieldTypes>
        static word fieldsName
        (
            const AlphaRhoFieldType& alphaRhoField,
            const AlphaRhoFieldTypes& ... alphaRhoFields
        );

        //- Return the name of the product of the fields associated with a
        //  source term (base condition for the above)
        template<class AlphaRhoFieldType>
        static word fieldsName(const AlphaRhoFieldType& alphaRhoField);


    // Constructors

        //- Construct from components
        LagrangianModel
        (
            const word& name,
            const LagrangianMesh& mesh
        );

        //- Disallow default bitwise copy construction
        LagrangianModel(const LagrangianModel&) = delete;

        //- Clone
        autoPtr<LagrangianModel> clone() const
        {
            NotImplemented;
            return autoPtr<LagrangianModel>(nullptr);
        }

        //- List construction class
        class iNew
        {
            const word& name_;

            const LagrangianMesh& mesh_;

        public:

            iNew(const word& name, const LagrangianMesh& mesh)
            :
                name_(name),
                mesh_(mesh)
            {}

            autoPtr<LagrangianModel> operator()(Istream& is) const
            {
                return autoPtr<LagrangianModel>
                (
                    LagrangianModel::New(name_, mesh_, dictionary(is))
                );
            }
        };


    //- Selector
    static autoPtr<LagrangianModel> New
    (
        const word& name,
        const LagrangianMesh& mesh,
        const dictionary& modelDict
    );


    //- Destructor
    virtual ~LagrangianModel();


    // Member Functions

        // Access

            //- The source name
            inline const word& name() const;

            //- The database
            inline const objectRegistry& db() const;

            //- The mesh
            inline const LagrangianMesh& mesh() const;


        // Checks

            //- Return the list of fields for which the LagrangianModel adds
            //  source term to the transport equation
            virtual wordList addSupFields() const;

            //- Return true if the LagrangianModel adds a source term to the
            //  given field's transport equation
            virtual bool addsSupToField(const word&) const;

            //- Return true if the LagrangianModels adds a source term to the
            //  given field's transport equation
            template<class Type, template<class> class PrimitiveField>
            bool addsSupToField
            (
                const LagrangianSubField<Type, PrimitiveField>& field
            ) const;


        //- Do post construction steps which require access to other models
        virtual void postConstruct();

        //- Correct the LagrangianModel
        virtual void correct();

        //- Identify elements in the Lagrangian mesh which are to be
        //  instantaneously modified or removed
        virtual void preModify
        (
            const LagrangianMesh& mesh,
            DynamicList<elementModification>& elementModifications
        ) const;

        //- Instantaneously modify and/or create and remove elements in the
        //  Lagrangian mesh
        virtual LagrangianSubMesh modify
        (
            LagrangianMesh& mesh,
            const LagrangianSubMesh& modifiedMesh
        ) const;

        //- Solve equations and/or update continually changing properties
        virtual void calculate
        (
            const LagrangianSubScalarField& deltaT,
            const bool final
        );


        // Sources

            //- Add a fractional source term
            virtual void addSup
            (
                const LagrangianSubScalarField& deltaT,
                LagrangianSubScalarField& S
            ) const;

            //- Add a source term to an equation
            FOR_ALL_FIELD_TYPES(DEFINE_LAGRANGIAN_MODEL_ADD_FIELD_SUP)

            //- Add a source term to a mass-weighted equation
            FOR_ALL_FIELD_TYPES(DEFINE_LAGRANGIAN_MODEL_ADD_M_FIELD_SUP)


        // Mesh changes

            //- Update topology using the given map
            virtual void topoChange(const polyTopoChangeMap&);

            //- Update from another mesh using the given map
            virtual void mapMesh(const polyMeshMap&);

            //- Redistribute or update using the given distribution map
            virtual void distribute(const polyDistributionMap&);


        // IO

            //- Read dictionary
            virtual bool read(const dictionary& modelDict);

            //- Write data
            virtual bool write(const bool write) const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const LagrangianModel&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "LagrangianModelI.H"

#ifdef NoRepository
    #include "LagrangianModelTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
