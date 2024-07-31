/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2024 OpenFOAM Foundation
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
    Foam::fv::homogeneousLiquidPhaseSeparation

Description
    Model for the homogeneous nucleation of a solid or liquid phase separating
    out of a liquid solution

Usage
    Example usage:
    \verbatim
    homogeneousLiquidPhaseSeparation
    {
        type            homogeneousLiquidPhaseSeparation;
        libs            ("libmultiphaseEulerFvModels.so");

        // Phases between which the transfer occurs. The first phase is the
        // solution, and the second is the precipitate.
        phases          (liquid sugar);

        // The specie that is condensing
        specie          C2H12O6;

        // Linearise the latent heat contribution into the energy equation?
        energySemiImplicit no;

        // Solubility given in mass of solute per mass of solvent
        solubility      constant 0.9;
    }
    \endverbatim

SourceFiles
    homogeneousLiquidPhaseSeparation.C

\*---------------------------------------------------------------------------*/

#ifndef homogeneousLiquidPhaseSeparation_H
#define homogeneousLiquidPhaseSeparation_H

#include "singleComponentPhaseChange.H"
#include "nucleation.H"
#include "dispersedPhaseInterface.H"
#include "saturationPressureModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
              Class homogeneousLiquidPhaseSeparation Declaration
\*---------------------------------------------------------------------------*/

class homogeneousLiquidPhaseSeparation
:
    public singleComponentPhaseChange,
    public nucleation
{
private:

    // Private Data

        //- Phase system
        const phaseSystem& fluid_;

        //- Interface
        const dispersedPhaseInterface interface_;

        //- Diameter of nucleated clusters
        volScalarField::Internal d_;

        //- Phase change rate, per unit volume of the solution
        volScalarField::Internal mDotByAlphaSolution_;

        //- Solubility curve
        autoPtr<Function1<scalar>> solubilityCurve_;


    // Private Member Functions

        //- Non-virtual read
        void readCoeffs(const dictionary& dict);

        //- Evaluate the solubility curve at the given temperature to return
        //  the saturated mass fraction of solute
        tmp<volScalarField::Internal> YSat
        (
            const volScalarField::Internal& T
        ) const;


public:

    //- Runtime type information
    TypeName("homogeneousLiquidPhaseSeparation");


    // Constructors

        //- Construct from explicit source name and mesh
        homogeneousLiquidPhaseSeparation
        (
            const word& name,
            const word& modelType,
            const fvMesh& mesh,
            const dictionary& dict
        );


    // Member Functions

        // Sources

            //- Return the diameter of nuclei
            virtual tmp<DimensionedField<scalar, volMesh>> d() const;

            //- Return the number rate at which nuclei are generated
            virtual tmp<DimensionedField<scalar, volMesh>> nDot() const;

            //- Return the mass transfer rate
            virtual tmp<DimensionedField<scalar, volMesh>> mDot() const;

            //- Override the compressible continuity equation to add
            //  linearisation w.r.t alpha
            void addSup
            (
                const volScalarField& alpha,
                const volScalarField& rho,
                fvMatrix<scalar>& eqn
            ) const;

            //- Override the compressible mass fraction equation to add
            //  linearisation w.r.t the mass fraction
            void addSup
            (
                const volScalarField& alpha,
                const volScalarField& rho,
                const volScalarField& Yi,
                fvMatrix<scalar>& eqn
            ) const;


        //- Correct the fvModel
        //  e.g. solve equations, update model, for film, Lagrangian etc.
        virtual void correct();


        // IO

            //- Read source dictionary
            virtual bool read(const dictionary& dict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
