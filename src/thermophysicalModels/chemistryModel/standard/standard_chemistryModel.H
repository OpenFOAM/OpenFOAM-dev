/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2022-2026 OpenFOAM Foundation
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
    Foam::chemistryModels::standard

Description
    Non-templated Base class for Foam::chemistryModels::Standard

    With run-time selected ODE solver and optional cellZone to limit the region
    in which the chemistry is integrated to speed-up large cases in which the
    chemistry is active in a small sub-region.

    Example of \c standard \c chemistryModel specification in
    constant/chemistryProperties using the \c seulex ODE solver to integrate
    the included reaction system:
    \verbatim
    type            standard;

    initialChemicalTimeStep 1e-7;

    ode
    {
        solver          seulex;
        absTol          1e-8;
        relTol          1e-1;
    }

    #include "reactions"
    \endverbatim

See also
    Foam::chemistryModel
    Foam::chemistryModels::Standard
    Foam::generatedCellZone

SourceFiles
    standard_chemistryModelI.H
    standard_chemistryModel.C

\*---------------------------------------------------------------------------*/

#ifndef standard_chemistryModel_H
#define standard_chemistryModel_H

#include "chemistryModel.H"
#include "ODESolver.H"
#include "specieThermosRef.H"
#include "reactionList.H"
#include "generatedCellZone.H"
#include "DynamicField.H"
#include "OFstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class chemistryReductionMethod;
class chemistryTabulationMethod;

namespace chemistryModels
{

/*---------------------------------------------------------------------------*\
                           Class standard Declaration
\*---------------------------------------------------------------------------*/

class standard
:
    public chemistryModel,
    public ODESystem
{
    // Private Classes

        //- Class to define scope of reaction evaluation. Runs pre-evaluate
        //  hook on all reactions on construction and post-evaluate on
        //  destruction.
        class reactionEvaluationScope
        {
            const standard& chemistry_;

        public:

            reactionEvaluationScope(const standard& chemistry)
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


    // Private Member Functions

        //- Return the species' molar weights in a list
        static scalarList specieWs(const fluidMulticomponentThermo& thermo);

        //- Return the species' formation enthalpies in a list
        static scalarList specieHfs(const fluidMulticomponentThermo& thermo);

        //- Solve the reaction system for the given time step
        //  of given type and return the characteristic time
        //  Variable number of species added
        template<class DeltaTType>
        scalar solve(const DeltaTType& deltaT);


protected:

    // Protected Data

        //- Enumeration for the type of Jacobian to be calculated by the
        //  chemistryModel
        enum class jacobianType
        {
            fast,
            exact
        };

        //- Jacobian type names
        static const NamedEnum<jacobianType, 2> jacobianTypeNames;

        //- Reference to the field of specie mass fractions
        const PtrList<volScalarField>& Yvf_;

        //- Number of species
        label nSpecie_;

        //- Reactions
        const reactionList reactions_;

        //- Reduction method
        autoPtr<chemistryReductionMethod> reductionPtr_;

        //- Reduction method reference
        chemistryReductionMethod& reduction_;

        //- Is chemistry reduction active?
        bool reduced_;

        //- Tabulation method
        autoPtr<chemistryTabulationMethod> tabulationPtr_;

        //- Tabulation method reference
        chemistryTabulationMethod& tabulation_;

        //- Optional cellZone to limit where the chemistry is integrated
        generatedCellZone zone_;

        //- Local cache of the species' molecular weights
        const scalarList specieWs_;

        //- Local cache of the species' formation enthalpies
        const scalarList specieHfs_;

        //- Temporary map from complete to simplified concentration fields
        //  c -> sc
        mutable List<label> c2s_;

        //- Temporary map from simplified to complete concentration fields
        //  sc -> c
        mutable DynamicList<label> sToc_;

        //- ODE solver
        mutable autoPtr<ODESolver> odeSolver_;

        //- ODE solver integration state variable
        mutable scalarField cT_;

        //- Constant pressure for the current ODE integration
        mutable scalar p_;

        //- List of reaction rate per specie [kg/m^3/s]
        PtrList<volInternalScalarField> RR_;

        //- Temporary mass fraction field
        mutable scalarField Y_;

        //- Temporary simplified mechanism mass fraction field
        mutable DynamicField<scalar> sY_;

        //- Temporary reaction state
        mutable reaction::state s_;

        //- Specie-temperature workspace fields
        mutable FixedList<scalarField, 3> YTWork_;

        //- Specie-temperature workspace matrices
        mutable FixedList<scalarSquareMatrix, 2> YTYTWork_;

        //- Log file for average time spent solving the chemistry
        autoPtr<OFstream> cpuSolveFile_;

        //- Switch to select performance logging
        Switch log_;

        //- Switch to enable per-cell CPU load caching for load-balancing
        Switch cpuLoad_;


public:

    //- Runtime type information
    TypeName("standard");


    // Constructors

        //- Construct from thermo
        standard
        (
            const fluidMulticomponentThermo& thermo,
            const specieThermosRef& specieThermos
        );

        //- Disallow default bitwise copy construction
        standard(const standard&) = delete;


    //- Destructor
    virtual ~standard();


    // Member Functions

        //- The reactions
        inline const PtrList<reaction>& reactions() const;

        //- Create and return a TDAC log file of the given name
        inline autoPtr<OFstream> logFile(const word& name) const;

        //- The number of species
        inline virtual label nSpecie() const;

        //- The number of reactions
        inline virtual label nReaction() const;

        //- Return the name of reactioni
        virtual inline const word& reactionName(const label reactioni) const;

        //- Allow the reduction method to reset the number of species
        inline void setNSpecie(const label newNs);

        //- Number of ODE's to solve
        inline virtual label nEqns() const;

        //- Return a reference to the list of mass fraction fields
        inline const PtrList<volScalarField>& Y() const;

        //- Return true if mechanism reduction is active
        inline bool reduced() const;

        //- Return the index in the complete set of species
        //  corresponding to the index si in the simplified set of species
        inline label sToc(const label si) const;

        //- Return the index in the simplified set of species
        //  corresponding to the index si in the complete set of species
        inline label c2s(const label ci) const;


        // Overrides to chemistryModel functions

            //- Return reaction rates of the species [kg/m^3/s]
            virtual inline const PtrList<volInternalScalarField>& RR() const;

            //- Return the rate of reactioni [kmol/m^3/s]
            virtual tmp<volInternalScalarField> reactionRR
            (
                const label reactioni
            ) const;

            //- Return reaction rates of the species in reactioni [kg/m^3/s]
            virtual PtrList<volInternalScalarField> specieReactionRR
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
            virtual tmp<volInternalScalarField> tc() const;

            //- Return the heat release rate [kg/m/s^3]
            virtual tmp<volInternalScalarField> Qdot() const;


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
            ) const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const standard&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace chemistryModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "standard_chemistryModelI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
