/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2021 OpenFOAM Foundation
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

#include "basicChemistryModel.H"
#include "ReactionList.H"
#include "ODESystem.H"
#include "volFields.H"
#include "multiComponentMixture.H"
#include "chemistryReductionMethod.H"
#include "chemistryTabulationMethod.H"
#include "DynamicField.H"
#include "OFstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class chemistryModel Declaration
\*---------------------------------------------------------------------------*/

template<class ThermoType>
class chemistryModel
:
    public basicChemistryModel,
    public ODESystem
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

        //- Type of the Jacobian to be calculated
        const jacobianType jacobianType_;

        //- Reference to the field of specie mass fractions
        const PtrList<volScalarField>& Yvf_;

        //- Reference to the multi component mixture
        const multiComponentMixture<ThermoType>& mixture_;

        //- Thermodynamic data of the species
        const PtrList<ThermoType>& specieThermos_;

        //- Reactions
        const ReactionList<ThermoType> reactions_;

        //- Number of species
        label nSpecie_;

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

        //- Temporary map from complete to simplified concentration fields
        //  c -> sc
        List<label> cTos_;

        //- Temporary map from simplified to complete concentration fields
        //  sc -> c
        DynamicList<label> sToc_;

        //- Mechanism reduction method
        autoPtr<chemistryReductionMethod<ThermoType>> mechRedPtr_;

        //- Mechanism reduction method reference
        chemistryReductionMethod<ThermoType>& mechRed_;

        bool mechRedActive_;

        //- Tabulation method
        autoPtr<chemistryTabulationMethod<ThermoType>> tabulationPtr_;

        //- Tabulation method reference
        chemistryTabulationMethod<ThermoType>& tabulation_;

        //- Log file for average time spent solving the chemistry
        autoPtr<OFstream> cpuSolveFile_;


    // Private Member Functions

        //- Write access to chemical source terms
        //  (e.g. for multi-chemistry model)
        inline PtrList<volScalarField::Internal>& RR();

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
        chemistryModel(const fluidReactionThermo& thermo);

        //- Disallow default bitwise copy construction
        chemistryModel(const chemistryModel&) = delete;


    //- Destructor
    virtual ~chemistryModel();


    // Member Functions

        //- Create and return a TDAC log file of the given name
        inline autoPtr<OFstream> logFile(const word& name) const;

        //- Return reference to the mixture
        inline const multiComponentMixture<ThermoType>& mixture() const;

        //- The reactions
        inline const PtrList<Reaction<ThermoType>>& reactions() const;

        //- Thermodynamic data of the species
        inline const PtrList<ThermoType>& specieThermos() const;

        //- The number of species
        virtual inline label nSpecie() const;

        //- The number of reactions
        virtual inline label nReaction() const;

        //- Calculates the reaction rates
        virtual void calculate();


        // Chemistry model functions (overriding abstract functions in
        // basicChemistryModel.H)

            //- Return const access to the chemical source terms for specie, i
            inline const volScalarField::Internal& RR
            (
                const label i
            ) const;

            //- Return non const access to chemical source terms [kg/m^3/s]
            virtual volScalarField::Internal& RR
            (
                const label i
            );

            //- Return reaction rate of the speciei in reactionI
            virtual tmp<volScalarField::Internal> calculateRR
            (
                const label reactionI,
                const label speciei
            ) const;

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


        // ODE functions (overriding abstract functions in ODE.H)

            //- Number of ODE's to solve
            inline virtual label nEqns() const;

            virtual void derivatives
            (
                const scalar t,
                const scalarField& YTp,
                const label li,
                scalarField& dYTpdt
            ) const;

            virtual void jacobian
            (
                const scalar t,
                const scalarField& YTp,
                const label li,
                scalarField& dYTpdt,
                scalarSquareMatrix& J
            ) const;

            virtual void solve
            (
                scalar& p,
                scalar& T,
                scalarField& Y,
                const label li,
                scalar& deltaT,
                scalar& subDeltaT
            ) const = 0;


        //- Return a reference to the list of mass fraction fields
        inline const PtrList<volScalarField>& Y() const;


        // Mechanism reduction access functions

            //- Return access to the mechanism reduction method
            //  Needed by ISAT
            inline const chemistryReductionMethod<ThermoType>& mechRed() const;

            //- Allow the reduction method to reset the number of species
            inline void setNSpecie(const label newNs);

            //- Return true if specie i is active
            inline bool active(const label i) const;

            //- Set specie i active
            inline void setActive(const label i);

            //- Return the index in the complete set of species
            //  corresponding to the index si in the simplified set of species
            inline label sToc(const label si) const;

            //- Return the index in the simplified set of species
            //  corresponding to the index si in the complete set of species
            inline label cTos(const label ci) const;


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
