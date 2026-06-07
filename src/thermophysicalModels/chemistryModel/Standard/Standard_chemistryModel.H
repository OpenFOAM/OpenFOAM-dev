/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2026 OpenFOAM Foundation
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
    Foam::chemistryModels::Standard

Description
    Extends Foam::chemistryModels::standard

    adding a thermo package, and ODE functions.

    Integrates a standard OpenFOAM reaction system and evaluation of chemical
    source terms with optional support for mechanism reduction and tabulation.

SourceFiles
    Standard_chemistryModelI.H
    standard_chemistryModel.C

\*---------------------------------------------------------------------------*/

#ifndef Standard_chemistryModel_H
#define Standard_chemistryModel_H

#include "standard_chemistryModel.H"
#include "ReactionList.H"
#include "ODESystem.H"
#include "volFields.H"
#include "multicomponentMixture.H"
#include "chemistryReductionMethod.H"
#include "chemistryTabulationMethod.H"
#include "DynamicField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace chemistryModels
{

/*---------------------------------------------------------------------------*\
                           Class Standard Declaration
\*---------------------------------------------------------------------------*/

template<class ThermoType>
class Standard
:
    public standard
{
    // Private classes

        //- Class to define scope of reaction evaluation. Runs pre-evaluate
        //  hook on all reactions on construction and post-evaluate on
        //  destruction.
        class reactionEvaluationScope
        {
            const Standard<ThermoType>& chemistry_;

        public:

            reactionEvaluationScope
            (
                const Standard<ThermoType>& chemistry
            )
            :
                chemistry_(chemistry)
            {
                forAll(chemistry_.reactions_, i)
                {
                    chemistry_.reactions_[i].preEvaluate();
                }
            }

            ~reactionEvaluationScope()
            {
                forAll(chemistry_.reactions_, i)
                {
                    chemistry_.reactions_[i].postEvaluate();
                }
            }
        };


    // Private data

        //- Switch to select performance logging
        Switch log_;

        //- Switch to enable per-cell CPU load caching for load-balancing
        Switch cpuLoad_;

        //- Type of the Jacobian to be calculated
        const jacobianType jacobianType_;

        //- Reference to the multi component mixture
        const multicomponentMixture<ThermoType>& mixture_;

        //- Thermodynamic data of the species
        const PtrList<ThermoType>& specieThermos_;

        //- Reactions
        const ReactionList<ThermoType> reactions_;

        //- List of reaction rate per specie [kg/m^3/s]
        PtrList<volScalarField::Internal> RR_;

        //- Temporary mass fraction field
        mutable scalarField Y_;

        //- Temporary simplified mechanism mass fraction field
        DynamicField<scalar> sY_;

        //- Temporary concentration field
        mutable scalarField c_;

        //- Temporary simplified mechanism concentration field
        DynamicField<scalar> sc_;

        //- Specie-temperature-pressure workspace fields
        mutable FixedList<scalarField, 5> YTpWork_;

        //- Specie-temperature-pressure workspace matrices
        mutable FixedList<scalarSquareMatrix, 2> YTpYTpWork_;

        //- Mechanism reduction method
        autoPtr<chemistryReductionMethod<ThermoType>> mechRedPtr_;

        //- Mechanism reduction method reference
        chemistryReductionMethod<ThermoType>& mechRed_;

        //- Tabulation method
        autoPtr<chemistryTabulationMethod> tabulationPtr_;

        //- Tabulation method reference
        chemistryTabulationMethod& tabulation_;

        //- Log file for average time spent solving the chemistry
        autoPtr<OFstream> cpuSolveFile_;


    // Private Member Functions

        //- Solve the reaction system for the given time step
        //  of given type and return the characteristic time
        //  Variable number of species added
        template<class DeltaTType>
        scalar solve(const DeltaTType& deltaT);


public:

    //- Runtime type information
    TypeName("standard");


    // Constructors

        //- Construct from thermo
        Standard(const fluidMulticomponentThermo& thermo);

        //- Disallow default bitwise copy construction
        Standard(const Standard&) = delete;


    //- Destructor
    virtual ~Standard();


    // Member Functions

        // Access

            //- Return reference to the mixture
            inline const multicomponentMixture<ThermoType>& mixture() const;

            //- The reactions
            inline const PtrList<Reaction<ThermoType>>& reactions() const;

            //- Thermodynamic data of the species
            inline const PtrList<ThermoType>& specieThermos() const;


        // Overrides to basicChemistryModel functions

            //- The number of reactions
            virtual inline label nReaction() const;

            //- Return reaction rates of the species [kg/m^3/s]
            virtual inline const PtrList<volScalarField::Internal>& RR() const;

            //- Return the name of reactioni
            virtual inline const word& reactionName
            (
                const label reactioni
            ) const;

            //- Return the rate of reactioni [kmol/m^3/s]
            virtual tmp<volScalarField::Internal> reactionRR
            (
                const label reactioni
            ) const;

            //- Return reaction rates of the species in reactioni [kg/m^3/s]
            virtual PtrList<volScalarField::Internal> specieReactionRR
            (
                const label reactioni
            ) const;

            //- Calculates the reaction rates
            virtual void calculate();

            //- Solve the reaction system for the given time step
            //  and return the characteristic time
            virtual scalar solve(const scalar deltaT);

            //- Solve the reaction system for the given time step
            //  and return the characteristic time
            virtual scalar solve(const scalarField& deltaT);

            //- Return the chemical time scale
            virtual tmp<volScalarField> tc() const;

            //- Return the heat release rate [kg/m/s^3]
            virtual tmp<volScalarField> Qdot() const;


        // Overrides to ODESystem

            //- Calculate the ODE derivatives
            virtual void derivatives
            (
                const scalar t,
                const scalarField& YTp,
                const label li,
                scalarField& dYTpdt
            ) const;

            //- Calculate the ODE jacobian
            virtual void jacobian
            (
                const scalar t,
                const scalarField& YTp,
                const label li,
                scalarField& dYTpdt,
                scalarSquareMatrix& J
            ) const;


        // ODE solution functions

            //- Solve the ODE system
            virtual void solve
            (
                scalar& p,
                scalar& T,
                scalarField& Y,
                const label li,
                scalar& deltaT,
                scalar& subDeltaT
            ) const = 0;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const Standard&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace chemistryModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "Standard_chemistryModelI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "Standard_chemistryModel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
