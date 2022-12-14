/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2022 OpenFOAM Foundation
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
    Foam::chemistryModel

Description
    Extends base chemistry model by adding a thermo package, and ODE functions.
    Introduces chemistry equation system and evaluation of chemical source terms
    with optional support for TDAC mechanism reduction and tabulation.

    References:
    \verbatim
        Contino, F., Jeanmart, H., Lucchini, T., & D’Errico, G. (2011).
        Coupling of in situ adaptive tabulation and dynamic adaptive chemistry:
        An effective method for solving combustion in engine simulations.
        Proceedings of the Combustion Institute, 33(2), 3057-3064.

        Contino, F., Lucchini, T., D'Errico, G., Duynslaegher, C.,
        Dias, V., & Jeanmart, H. (2012).
        Simulations of advanced combustion modes using detailed chemistry
        combined with tabulation and mechanism reduction techniques.
        SAE International Journal of Engines,
        5(2012-01-0145), 185-196.

        Contino, F., Foucher, F., Dagaut, P., Lucchini, T., D’Errico, G., &
        Mounaïm-Rousselle, C. (2013).
        Experimental and numerical analysis of nitric oxide effect on the
        ignition of iso-octane in a single cylinder HCCI engine.
        Combustion and Flame, 160(8), 1476-1483.

        Contino, F., Masurier, J. B., Foucher, F., Lucchini, T., D’Errico, G., &
        Dagaut, P. (2014).
        CFD simulations using the TDAC method to model iso-octane combustion
        for a large range of ozone seeding and temperature conditions
        in a single cylinder HCCI engine.
        Fuel, 137, 179-184.
    \endverbatim

SourceFiles
    chemistryModelI.H
    chemistryModel.C

\*---------------------------------------------------------------------------*/

#ifndef chemistryModel_H
#define chemistryModel_H

#include "odeChemistryModel.H"
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

/*---------------------------------------------------------------------------*\
                     Class chemistryModel Declaration
\*---------------------------------------------------------------------------*/

template<class ThermoType>
class chemistryModel
:
    public odeChemistryModel
{
    // Private classes

        //- Class to define scope of reaction evaluation. Runs pre-evaluate
        //  hook on all reactions on construction and post-evaluate on
        //  destruction.
        class reactionEvaluationScope
        {
            const chemistryModel<ThermoType>& chemistry_;

        public:

            reactionEvaluationScope
            (
                const chemistryModel<ThermoType>& chemistry
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

        //- Switch to enable loadBalancing performance logging
        Switch loadBalancing_;

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
    TypeName("chemistryModel");


    // Constructors

        //- Construct from thermo
        chemistryModel(const fluidMulticomponentThermo& thermo);

        //- Disallow default bitwise copy construction
        chemistryModel(const chemistryModel&) = delete;


    //- Destructor
    virtual ~chemistryModel();


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

            //- Return reaction rates of the species in reactioni [kg/m^3/s]
            virtual PtrList<volScalarField::Internal> reactionRR
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


        // Mechanism reduction functions

            //- Return true if specie i is active
            inline bool active(const label i) const;

            //- Set specie i active
            inline void setActive(const label i);


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const chemistryModel&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "chemistryModelI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "chemistryModel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
