/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2021 OpenFOAM Foundation
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
    Foam::fvModel

Description
    Finite volume model abstract base class.

SourceFiles
    fvModel.C

\*---------------------------------------------------------------------------*/

#ifndef fvModel_H
#define fvModel_H

#include "fvMatricesFwd.H"
#include "volFieldsFwd.H"
#include "dictionary.H"
#include "dimensionSet.H"
#include "fvModelM.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class fvMesh;
class mapPolyMesh;

/*---------------------------------------------------------------------------*\
                           Class fvModel Declaration
\*---------------------------------------------------------------------------*/

class fvModel
{
    // Private Member Data

        //- Model name
        const word name_;

        //- Model type
        const word modelType_;

        //- Reference to the mesh database
        const fvMesh& mesh_;

        //- Top level source dictionary
        dictionary dict_;

        //- Dictionary containing source coefficients
        dictionary coeffs_;


protected:

    // Protected Member Functions

        //- Add a source term to an equation
        template<class Type>
        void addSupType
        (
            fvMatrix<Type>& eqn,
            const word& fieldName
        ) const;

        //- Add a source term to a compressible equation
        template<class Type>
        void addSupType
        (
            const volScalarField& rho,
            fvMatrix<Type>& eqn,
            const word& fieldName
        ) const;

        //- Add a source term to a phase equation
        template<class Type>
        void addSupType
        (
            const volScalarField& alpha,
            const volScalarField& rho,
            fvMatrix<Type>& eqn,
            const word& fieldName
        ) const;


        //- Return source for equation with specified name and dimensions
        template<class Type, class ... AlphaRhoFieldTypes>
        tmp<fvMatrix<Type>> source
        (
            const GeometricField<Type, fvPatchField, volMesh>& field,
            const word& fieldName,
            const dimensionSet& ds,
            const AlphaRhoFieldTypes& ... alphaRhos
        ) const;


public:

    //- Runtime type information
    TypeName("fvModel");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            fvModel,
            dictionary,
            (
                const word& name,
                const word& modelType,
                const dictionary& dict,
                const fvMesh& mesh
            ),
            (name, modelType, dict, mesh)
        );


    // Static Member Functions

        //- Return the dimensions of the matrix of a source term
        template
        <
            class Type,
            class AlphaRhoFieldType,
            class ... AlphaRhoFieldTypes
        >
        inline static dimensionSet sourceDims
        (
            const GeometricField<Type, fvPatchField, volMesh>& field,
            const dimensionSet& ds,
            const AlphaRhoFieldType& alphaRho,
            const AlphaRhoFieldTypes& ... alphaRhos
        );

        //- Return the dimensions of the matrix of a source term (base
        //  condition for the above)
        template<class Type>
        inline static dimensionSet sourceDims
        (
            const GeometricField<Type, fvPatchField, volMesh>& field,
            const dimensionSet& ds
        );


    // Constructors

        //- Construct from components
        fvModel
        (
            const word& name,
            const word& modelType,
            const dictionary& dict,
            const fvMesh& mesh
        );

        //- Return clone
        autoPtr<fvModel> clone() const
        {
            NotImplemented;
            return autoPtr<fvModel>(nullptr);
        }

        //- Return pointer to new fvModel object created
        //  on the freestore from an Istream
        class iNew
        {
            //- Reference to the mesh
            const fvMesh& mesh_;

            const word& name_;

        public:

            iNew
            (
                const fvMesh& mesh,
                const word& name
            )
            :
                mesh_(mesh),
                name_(name)
            {}

            autoPtr<fvModel> operator()(Istream& is) const
            {
                // const word name(is);
                const dictionary dict(is);

                return autoPtr<fvModel>
                (
                    fvModel::New(name_, dict, mesh_)
                );
            }
        };


    // Selectors

        //- Return a reference to the selected fvModel
        static autoPtr<fvModel> New
        (
            const word& name,
            const dictionary& dict,
            const fvMesh& mesh
        );


    //- Destructor
    virtual ~fvModel();


    // Member Functions

        // Access

            //- Return const access to the source name
            inline const word& name() const;

            //- Return const access to the mesh database
            inline const fvMesh& mesh() const;

            //- Return dictionary
            inline const dictionary& coeffs() const;


        // Checks

            //- Return the list of fields for which the fvModel adds source term
            //  to the transport equation
            virtual wordList addSupFields() const;

            //- Return true if the fvModel adds a source term to the given
            //  field's transport equation
            virtual bool addsSupToField(const word& fieldName) const;


        // Sources

            //- Add a source term to an equation
            FOR_ALL_FIELD_TYPES(DEFINE_FV_MODEL_ADD_SUP);

            //- Add a source term to a compressible equation
            FOR_ALL_FIELD_TYPES(DEFINE_FV_MODEL_ADD_RHO_SUP);

            //- Add a source term to a phase equation
            FOR_ALL_FIELD_TYPES(DEFINE_FV_MODEL_ADD_ALPHA_RHO_SUP);

            //- Return source for an equation
            template<class Type>
            tmp<fvMatrix<Type>> source
            (
                const GeometricField<Type, fvPatchField, volMesh>& field
            ) const;

            //- Return source for an equation with a specified name
            template<class Type>
            tmp<fvMatrix<Type>> source
            (
                const GeometricField<Type, fvPatchField, volMesh>& field,
                const word& fieldName
            ) const;

            //- Return source for a compressible equation
            template<class Type>
            tmp<fvMatrix<Type>> source
            (
                const volScalarField& rho,
                const GeometricField<Type, fvPatchField, volMesh>& field
            ) const;

            //- Return source for a compressible equation with a specified name
            template<class Type>
            tmp<fvMatrix<Type>> source
            (
                const volScalarField& rho,
                const GeometricField<Type, fvPatchField, volMesh>& field,
                const word& fieldName
            ) const;

            //- Return source for a phase equation
            template<class Type>
            tmp<fvMatrix<Type>> source
            (
                const volScalarField& alpha,
                const volScalarField& rho,
                const GeometricField<Type, fvPatchField, volMesh>& field
            ) const;

            //- Return source for a phase equation with a specified name
            template<class Type>
            tmp<fvMatrix<Type>> source
            (
                const volScalarField& alpha,
                const volScalarField& rho,
                const GeometricField<Type, fvPatchField, volMesh>& field,
                const word& fieldName
            ) const;


        // Mesh changes

            //- Prepare for mesh update
            virtual void preUpdateMesh();

            //- Update for mesh changes
            virtual void updateMesh(const mapPolyMesh&);

            //- Update for mesh motion
            virtual bool movePoints();


        //- Correct the fvModel
        //  e.g. solve equations, update model, for film, Lagrangian etc.
        virtual void correct();


        // IO

            //- Read source dictionary
            virtual bool read(const dictionary& dict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "fvModelI.H"

#ifdef NoRepository
    #include "fvModelTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
