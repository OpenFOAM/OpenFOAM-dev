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
    Foam::fv::cloud

Description
    Finite volume model that solves for the evolution of a cloud. Provides
    two-way coupling with a finite-volume carrier phase.

SourceFiles
    cloud_fvModel.C
    cloud_fvModelTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef cloud_fvModel_H
#define cloud_fvModel_H

#include "fvModel.H"
#include "coupled.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
                            Class cloud Declaration
\*---------------------------------------------------------------------------*/

class cloud
:
    public fvModel
{
    // Private Data

        //- The cloud
        autoPtr<Foam::cloud> cloudPtr_;

        //- Cloud reference (for convenience)
        const Foam::cloud& cloud_;

        //- Coupled cloud reference (for convenience)
        const clouds::coupled& coupledCloud_;


    // Private Member Functions

        // Sources

            //- Add a source term to an equation
            template<class Type>
            void addSupType
            (
                const VolField<Type>& field,
                fvMatrix<Type>& eqn
            ) const;

            //- Add a source term to a compressible equation
            template<class Type>
            void addSupType
            (
                const volScalarField& rho,
                const VolField<Type>& field,
                fvMatrix<Type>& eqn
            ) const;

            //- Add a source term to a phase equation
            template<class Type>
            void addSupType
            (
                const volScalarField& alpha,
                const volScalarField& rho,
                const VolField<Type>& field,
                fvMatrix<Type>& eqn
            ) const;


protected:

    // Protected Classes

        //- Empty template class used to template the constructor on the type
        //  of the cloud
        template<class Type>
        class Cloud
        {};


    // Protected Constructors

        //- Construct from components
        template<class Type>
        cloud
        (
            const word& name,
            const word& modelType,
            const fvMesh& mesh,
            const dictionary& dict,
            const Cloud<Type>&
        );


public:

    // Public Static Data

        //- Run-time type information
        TypeName("cloud");


    //- Destructor
    virtual ~cloud();


    // Member Functions

        // Checks

            //- Return true if the fvModel adds a source term to the given
            //  field's transport equation
            virtual bool addsSupToField(const word& fieldName) const;


        // Sources

            //- Add a source term to a field-less proxy equation
            virtual void addSup(fvMatrix<scalar>& eqn) const;

            //- Add a source term to an equation
            FOR_ALL_FIELD_TYPES(DEFINE_FV_MODEL_ADD_FIELD_SUP)

            //- Add a source term to a compressible equation
            FOR_ALL_FIELD_TYPES(DEFINE_FV_MODEL_ADD_RHO_FIELD_SUP)

            //- Add a source term to a phase equation
            FOR_ALL_FIELD_TYPES(DEFINE_FV_MODEL_ADD_ALPHA_RHO_FIELD_SUP)


        //- Solve the evolution of the cloud over the current time-step
        virtual void correct();


        // Mesh changes

            //- Prepare for mesh update
            virtual void preUpdateMesh();

            //- Update for mesh motion
            virtual bool movePoints();

            //- Update topology using the given map
            virtual void topoChange(const polyTopoChangeMap&);

            //- Update from another mesh using the given map
            virtual void mapMesh(const polyMeshMap&);

            //- Redistribute or update using the given distribution map
            virtual void distribute(const polyDistributionMap&);
};


/*---------------------------------------------------------------------------*\
                            Class Cloud Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class Cloud
:
    public cloud
{
public:

    // Public Static Data

        //- Run-time type information
        TypeName("Cloud");


    // Constructors

        //- Construct from components
        Cloud
        (
            const word& name,
            const word& modelType,
            const fvMesh& mesh,
            const dictionary& dict
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "cloud_fvModelTemplates.C"
#endif


#define makeCloudFvModel(Type)                                              \
                                                                            \
    typedef Cloud<clouds::Type> Type##CloudFvModel;                         \
                                                                            \
    defineTemplateTypeNameAndDebugWithName                                  \
    (                                                                       \
        Type##CloudFvModel,                                                 \
        (clouds::Type::typeName + "Cloud").c_str(),                         \
        0                                                                   \
    );                                                                      \
                                                                            \
    addToRunTimeSelectionTable                                              \
    (                                                                       \
        fvModel,                                                            \
        Type##CloudFvModel,                                                 \
        dictionary                                                          \
    );


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
