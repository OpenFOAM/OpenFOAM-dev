/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2017-2021 OpenFOAM Foundation
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
    Foam::diameterModels::populationBalanceModel

Description
    Class that solves the univariate population balance equation by means of
    a class method (also called sectional or discrete method). The internal
    coordinate is set to the particle volume, so the equation is based on
    a transport equation of the volume-based number density function. The
    discretisation is done using the fixed pivot technique of Kumar and
    Ramkrishna (1996). The source terms are written in a way that particle
    number and mass are preserved. Coalescence (coagulation), breakup, drift
    (growth and surface loss) as well as nucleation are supported.
    For the discrete breakup term two recipes are available, depending on the
    model choice. For models which state a total breakup rate and a separate
    daughter size distribution function, the formulation of Kumar and Ramkrishna
    (1996) is applied which is applicable for binary and multiple breakup
    events. The second formulation is given by Liao et al. (2018). It is useful
    for binary breakup models which give the breakup rate between a sizeGroup
    pair directly, without an explicit expression for the daughter size
    distribution. The drift term is implemented using a finite difference upwind
    scheme. Although it is diffusive, it ensures a stable and
    number-conservative solution.

    The implementation allows to split the population balance over multiple
    velocity fields using the capability of multiphaseEulerFoam to solve
    for n momentum equations. It is also possible to define multiple population
    balances, e.g. bubbles and droplets simultaneously.

    References:
    \verbatim
        Coalescence and breakup term formulation:
        Kumar, S., & Ramkrishna, D. (1996).
        On the solution of population balance equations by discretization-I. A
        fixed pivot technique.
        Chemical Engineering Science, 51(8), 1311-1332.
    \endverbatim

    \verbatim
        Binary breakup term formulation:
        Liao, Y., Oertel, R., Kriebitzsch, S., Schlegel, F., & Lucas, D. (2018).
        A discrete population balance equation for binary breakage.
        International Journal for Numerical Methods in Fluids, 87(4), 202-215.
    \endverbatim

Usage
    Example excerpt from a phaseProperties dictionary.
    \verbatim
    type populationBalanceTwoPhaseSystem;

    phases (air water);

    populationBalances (bubbles);

    air
    {
        type            purePhaseModel;
        diameterModel   velocityGroup;
        velocityGroupCoeffs
        {
            populationBalance    bubbles;

            shapeModel           constant;

            sizeGroups
            (
                f0{dSph  1.00e-3; value 0;}
                f1{dSph  1.08e-3; value 0;}
                f2{dSph  1.16e-3; value 0.25;}
                f3{dSph  1.25e-3; value 0.5;}
                f4{dSph  1.36e-3; value 0.25;}
                f5{dSph  1.46e-3; value 0;}
                ...
            );
        }

        residualAlpha   1e-6;
    }

    populationBalanceCoeffs
    {
        bubbles
        {
            continuousPhase water;

            coalescenceModels
            (
                hydrodynamic
                {
                    C 0.25;
                }
            );

            binaryBreakupModels
            ();

            breakupModels
            (
                exponential
                {
                    C 0.5;
                    exponent 0.01;
                    daughterSizeDistributionModel uniform;
                }
            );

            driftModels
            (
                densityChange{}
            );

            nucleationModels
            ();
        }
    }
    \endverbatim

See also
    Foam::diameterModels::sizeGroup
    Foam::diameterModels::velocityGroup

SourceFiles
    populationBalanceModel.C

\*---------------------------------------------------------------------------*/

#ifndef populationBalanceModel_H
#define populationBalanceModel_H

#include "sizeGroup.H"
#include "phasePair.H"
#include "pimpleControl.H"
#include "phaseDynamicMomentumTransportModel.H"
#include "HashPtrTable.H"
#include "Pair.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class phaseSystem;

namespace diameterModels
{

class coalescenceModel;
class breakupModel;
class binaryBreakupModel;
class driftModel;
class nucleationModel;

/*---------------------------------------------------------------------------*\
                   Class populationBalanceModel Declaration
\*---------------------------------------------------------------------------*/

class populationBalanceModel
:
    public regIOobject
{
private:

    // Private Typedefs

        typedef
            HashTable<autoPtr<phasePair>, phasePairKey, phasePairKey::hash>
            phasePairTable;

        typedef
            HashPtrTable<volScalarField, phasePairKey, phasePairKey::hash>
            pDmdtTable;


    // Private Data

        //- Reference to the phaseSystem
        const phaseSystem& fluid_;

        //- Interfacial mass transfer rates
        pDmdtTable& pDmdt_;

        //- Reference to the mesh
        const fvMesh& mesh_;

        //- Name of the populationBalance
        word name_;

        //- Dictionary
        dictionary dict_;

        //- Reference to pimpleControl
        const pimpleControl& pimple_;

        //- Continuous phase
        const phaseModel& continuousPhase_;

        //- velocityGroups belonging to this populationBalance
        UPtrList<velocityGroup> velocityGroups_;

        //- sizeGroups belonging to this populationBalance
        UPtrList<sizeGroup> sizeGroups_;

        //- List of unordered phasePairs in this populationBalance
        phasePairTable phasePairs_;

        //- sizeGroup boundaries
        PtrList<dimensionedScalar> v_;

        //- Section width required for binary breakup formulation
        PtrList<PtrList<dimensionedScalar>> delta_;

        //- Explicitly treated sources
        PtrList<volScalarField> Su_;

        //- Implicitly treated sources
        PtrList<volScalarField> Sp_;

        //- Field for caching sources
        volScalarField Sui_;

        //- Coalescence models
        PtrList<coalescenceModel> coalescenceModels_;

        //- Coalescence rate
        autoPtr<volScalarField> coalescenceRate_;

        //- Coalescence relevant size group pairs
        List<labelPair> coalescencePairs_;

        //- Breakup models
        PtrList<breakupModel> breakupModels_;

        //- Breakup rate
        autoPtr<volScalarField> breakupRate_;

        //- Binary breakup models
        PtrList<binaryBreakupModel> binaryBreakupModels_;

        //- Binary breakup rate
        autoPtr<volScalarField> binaryBreakupRate_;

        //- Binary breakup relevant size group pairs
        List<labelPair> binaryBreakupPairs_;

        //- Drift models
        PtrList<driftModel> drift_;

        //- Drift rate
        autoPtr<volScalarField> driftRate_;

        //- Zeroeth order models
        PtrList<nucleationModel> nucleation_;

        //- Zeroeth order rate
        autoPtr<volScalarField> nucleationRate_;

        //- Total void fraction
        autoPtr<volScalarField> alphas_;

        //- Mean Sauter diameter
        autoPtr<volScalarField> dsm_;

        //- Average velocity
        autoPtr<volVectorField> U_;

        //- Counter for interval between source term updates
        label sourceUpdateCounter_;


    // Private Member Functions

        void registerVelocityGroups();

        void registerSizeGroups(sizeGroup& group);

        void createPhasePairs();

        void precompute();

        void birthByCoalescence(const label j, const label k);

        void deathByCoalescence(const label i, const label j);

        void birthByBreakup(const label k, const label model);

        void deathByBreakup(const label i);

        void calcDeltas();

        void birthByBinaryBreakup(const label i, const label j);

        void deathByBinaryBreakup(const label j, const label i);

        void drift(const label i, driftModel& model);

        void nucleation(const label i, nucleationModel& model);

        void sources();

        void dmdt();

        void calcAlphas();

        tmp<volScalarField> calcDsm();

        void calcVelocity();

        //- Return whether the sources should be updated on this iteration
        bool updateSources();

        //- Return the number of corrections
        inline label nCorr() const;

        //- Return the interval at which the sources are updated
        inline label sourceUpdateInterval() const;

public:

    //- Runtime type information
    TypeName("populationBalanceModel");


    // Constructor

        populationBalanceModel
        (
            const phaseSystem& fluid,
            const word& name,
            HashPtrTable
            <
                volScalarField,
                phasePairKey,
                phasePairKey::hash
            >& pDmdt
        );

        //- Return clone
        autoPtr<populationBalanceModel> clone() const;

        //- Return a pointer to a new populationBalanceModel object created on
        //  freestore from Istream
        class iNew
        {
            const phaseSystem& fluid_;

            HashPtrTable<volScalarField, phasePairKey, phasePairKey::hash>&
                pDmdt_;

        public:

            iNew
            (
                const phaseSystem& fluid,
                HashPtrTable<volScalarField, phasePairKey, phasePairKey::hash>&
                    pDmdt
            )
            :
                fluid_(fluid),
                pDmdt_(pDmdt)
            {}

            autoPtr<populationBalanceModel> operator()(Istream& is) const
            {
                return autoPtr<populationBalanceModel>
                (
                    new populationBalanceModel(fluid_, word(is), pDmdt_)
                );
            }
        };


    //- Destructor
    virtual ~populationBalanceModel();

    // Member Functions

        //- Dummy write for regIOobject
        bool writeData(Ostream&) const;

        //- Return reference to the phaseSystem
        inline const phaseSystem& fluid() const;

        //- Return reference to the mesh
        inline const fvMesh& mesh() const;

        //- Return populationBalanceCoeffs dictionary
        inline const dictionary& dict() const;

        //- Return continuous phase
        inline const phaseModel& continuousPhase() const;

        //- Return the velocityGroups belonging to this populationBalance
        inline const UPtrList<velocityGroup>& velocityGroups() const;

        //- Return the sizeGroups belonging to this populationBalance
        inline const UPtrList<sizeGroup>& sizeGroups() const;

        //- Return list of unordered phasePairs in this populationBalance
        inline const phasePairTable& phasePairs() const;

        //- Return implicit source terms
        inline const volScalarField& Sp(const label i) const;

        //- Return coalescence relevant size group pairs
        inline const List<labelPair>& coalescencePairs() const;

        //- Return binary breakup relevant size group pairs
        inline const List<labelPair>& binaryBreakupPairs() const;

        //- Return total void of phases belonging to this populationBalance
        inline const volScalarField& alphas() const;

        //- Return average velocity
        inline const volVectorField& U() const;

        //- Returns true if both phases are velocity groups and
        //  belong to this populationBalance
        inline bool isVelocityGroupPair(const phasePair& pair) const;

        //- Return allocation coefficient
        const dimensionedScalar eta
        (
            const label i,
            const dimensionedScalar& v
        ) const;

        //- Return the surface tension coefficient between a given dispersed
        //  and the continuous phase
        const tmp<volScalarField> sigmaWithContinuousPhase
        (
            const phaseModel& dispersedPhase
        ) const;

        //- Return reference to turbulence model of the continuous phase
        const phaseCompressible::momentumTransportModel&
            continuousTurbulence() const;

        //- Solve the population balance equation
        void solve();

        //- Correct derived quantities
        void correct();
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace diameterModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "populationBalanceModelI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
