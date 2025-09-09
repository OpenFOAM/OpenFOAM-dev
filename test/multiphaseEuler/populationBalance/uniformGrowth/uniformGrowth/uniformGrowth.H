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
    Foam::fv::uniformGrowth

Description
    This fvModel applies a mass source which uniformly grows all the particles
    in a population balance by the same amount. This is not a physical model.
    It is designed for use with unit tests.

Usage
    Example usage:
    \verbatim
    airSource
    {
        type            uniformGrowth;
        libs            ("libuniformGrowth.so");

        populationBalance bubbles;

        massFlowRate    0.5;
    }
    \endverbatim

SourceFiles
    uniformGrowth.C

\*---------------------------------------------------------------------------*/

#ifndef uniformGrowth_H
#define uniformGrowth_H

#include "fvSpecificSource.H"
#include "populationBalanceModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
                        Class uniformGrowth Declaration
\*---------------------------------------------------------------------------*/

class uniformGrowth
:
    public fvSpecificSource
{
private:

    // Private Data

        //- Population balance model
        const diameterModels::populationBalanceModel& popBal_;

        //- Names of the phases
        const wordList phaseNames_;

        //- Names of the volume-fraction fields
        const wordList alphaNames_;

        //- Names of the density fields
        const wordList rhoNames_;

        //- Mass flow rate
        autoPtr<Function1<scalar>> massFlowRate_;


    // Private Member Functions

        //- Extract the phase names from the population balance model
        wordList phaseNames() const;

        //- Extract the volume-fraction names from the population balance model
        wordList alphaNames() const;

        //- Extract the density names from the population balance model
        wordList rhoNames() const;

        //- Non-virtual read
        void readCoeffs(const dictionary& dict);


        // Sources

            //- Add a source term to an equation
            template<class Type>
            void addSupType
            (
                const VolField<Type>& field,
                fvMatrix<Type>& eqn
            ) const;

            //- Add a source term to a scalar equation
            void addSupType
            (
                const volScalarField& rhoOrField,
                fvMatrix<scalar>& eqn
            ) const;

            //- Add a source term to a compressible equation
            template<class Type>
            void addSupType
            (
                const volScalarField& rho,
                const VolField<Type>& field,
                fvMatrix<Type>& eqn
            ) const;

            //- Add a source term to a compressible scalar equation
            void addSupType
            (
                const volScalarField& alphaOrRho,
                const volScalarField& rhoOrField,
                fvMatrix<scalar>& eqn
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


public:

    //- Runtime type information
    TypeName("uniformGrowth");


    // Constructors

        //- Construct from explicit source name and mesh
        uniformGrowth
        (
            const word& name,
            const word& modelType,
            const fvMesh& mesh,
            const dictionary& dict
        );


    // Member Functions

        // Sources

            //- Return true if the fvModel adds a source term to the given
            //  field's transport equation
            bool addsSupToField(const word& fieldName) const;

            //- Return the source value
            virtual tmp<DimensionedField<scalar, volMesh>> S
            (
                const word& fieldName
            ) const;

            //- Add a source term to a field-less proxy equation
            virtual void addSup(fvMatrix<scalar>& eqn) const;

            //- Add a source term to an equation
            FOR_ALL_FIELD_TYPES(DEFINE_FV_MODEL_ADD_FIELD_SUP);

            //- Add a source term to a compressible equation
            FOR_ALL_FIELD_TYPES(DEFINE_FV_MODEL_ADD_RHO_FIELD_SUP);

            //- Add a source term to a phase equation
            FOR_ALL_FIELD_TYPES(DEFINE_FV_MODEL_ADD_ALPHA_RHO_FIELD_SUP);


        //- Read source dictionary
        virtual bool read(const dictionary& dict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
