/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2017-2023 OpenFOAM Foundation
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
    Model for tracking the evolution of a dispersed phase size distribution due
    to coalescence (synonymous with coagulation, aggregation, agglomeration)
    and breakup events as well as density or phase changes. Provides an
    approximate solution of the population balance equation by means of a class
    method. The underlying theory is described in the article of Lehnigk et al.
    (2021).

    The size distribution, expressed through a volume-based number density
    function, is discretised using the fixot pivot technique of Kumar and
    Ramkrishna (1996). Thereby, the population balance equation is transformed
    into a series of transport equations for the particle (bubble, droplet)
    number concentrations in separate size classes that are coupled through
    their source terms. The discretisation is based on representative particle
    volumes, which are provided by the user through the corresponding sphere
    equivalent diameters.

    Since the representative volumes are fixed a priori and the total dispersed
    phase volume already available from solving the phase continuity equation,
    the model only determines the evolution of the individual size class
    fractions

    \f[
        f_{i,\varphi} = \frac{\alpha_{i,\varphi}}{\alpha_{\varphi}}\,,
    \f]

    where \f$\alpha_{i,\varphi}\f$ is the volume fraction of the size class and
    \f$\alpha_{\varphi}\f$ the total phase fraction of phase \f$\varphi\f$.

    The source terms are formulated such that the first and second moment of
    the distribution, i.e. the total particle number and volume, are conserved
    irrespective of the discretisation of the size domain. The treatment of
    particle breakup depends on the selected breakup submodels. For models
    which provide a total breakup frequency and a separate daughter size
    distribution function, the formulation provided Kumar and Ramkrishna (1996)
    is utilised, which is applicable both for binary and multiple breakup
    events. Currently, only field-independent daughter size distribution models
    are allowed. In case of binary breakup models that provide the breakup
    frequency between a size class pair directly, the formulation of Liao et
    al. (2018) is adopted, which is computationally more efficient compared to
    first extracting the field-dependent daughter size distribution and then
    consuming it in the formulation of Kumar and Ramkrishna. The source terms
    describing a drift of the size distribution through particle growth or
    shrinkage are derived using upwind differencing, thus ensuring conservation
    of the total particle number and volume. Note that due to the volume-based
    implementation, both density as well as phase change lead to a drift of the
    size distribution function. Further, users can specify multiple submodels
    for each mechanism, whose contributions are added up.

    The model also allows to distribute the size classes over multiple
    representative phases with identical physical properties that collectively
    define the dispersed phase. Thereby, size class fields can be transported
    with different velocity fields in order to account for the size dependency
    of the particle motion. A possible mass transfer between representative
    phases by means of coalescence, breakup and drift is taken into account.
    Similarly, the spatial evolution of secondary particle properties such as
    the particle surface area can be tracked.

    The key variable during a simulation is the Sauter diameter, which is
    computed from the size class fractions of the corresponding phase. The
    underlying size distribution can be extracted from the simulation using the
    functionObject 'sizeDistribution'. Integral and mean properties of a size
    distribution can be computed with the functionObject 'moments'.

    Verification cases for the population balance modeling functionality are
    provided in test/multiphase/multiphaseEulerFoam/populationBalance.

    References:
    \verbatim
        Lehnigk, R., Bainbridge, W., Liao, Y., Lucas, D., Niemi, T.,
        Peltola, J., & Schlegel, F. (2021).
        An open‐source population balance modeling framework for the simulation
        of polydisperse multiphase flows.
        AIChE Journal, 68(3), e17539.
    \endverbatim

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
    Excerpt from an exemplary phaseProperties dictionary:
    \verbatim
    type populationBalanceMultiphaseSystem;

    ...

    populationBalances (bubbles);

    air
    {
        type            pureIsothermalPhaseModel;

        diameterModel   velocityGroup;

        velocityGroupCoeffs
        {
            populationBalance    bubbles;

            shapeModel           spherical;

            sizeGroups
            (
                f1 {dSph  1e-3; value 1.0;}
                f2 {dSph  2e-3; value 0.0;}
                f3 {dSph  3e-3; value 0.0;}
                f4 {dSph  4e-3; value 0.0;}
                f5 {dSph  5e-3; value 0.0;}
                ...
            );
        }

        residualAlpha   1e-6;
    }

    ...

    populationBalanceCoeffs
    {
        bubbles
        {
            continuousPhase water;

            coalescenceModels
            (
                LehrMilliesMewes{}
            );

            binaryBreakupModels
            (
                LehrMilliesMewes{}
            );

            breakupModels
            ();

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
    Foam::PopulationBalancePhaseSystem
    Foam::diameterModels::sizeGroup
    Foam::diameterModels::velocityGroup
    Foam::diameterModels::SecondaryPropertyModel
    Foam::diameterModels::coalescenceModel
    Foam::diameterModels::breakupModel
    Foam::diameterModels::daughterSizeDistributionModel
    Foam::diameterModels::binaryBreakupModel
    Foam::diameterModels::driftModel
    Foam::diameterModels::nucleationModel
    Foam::functionObjects::sizeDistribution
    Foam::functionObjects::moments

SourceFiles
    populationBalanceModel.C

\*---------------------------------------------------------------------------*/

#ifndef populationBalanceModel_H
#define populationBalanceModel_H

#include "sizeGroup.H"
#include "phaseSystem.H"
#include "phaseCompressibleMomentumTransportModel.H"
#include "HashPtrTable.H"
#include "Pair.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
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
    // Private Data

        //- Reference to the phaseSystem
        const phaseSystem& fluid_;

        //- Interfacial mass transfer rates
        phaseSystem::dmdtfTable dmdtfs_;

        //- Reference to the mesh
        const fvMesh& mesh_;

        //- Name of the populationBalance
        word name_;

        //- Dictionary
        dictionary dict_;

        //- Continuous phase
        const phaseModel& continuousPhase_;

        //- Velocity groups belonging to this populationBalance
        HashTable<const velocityGroup*> velocityGroupPtrs_;

        //- Size groups belonging to this populationBalance
        UPtrList<sizeGroup> sizeGroups_;

        //- sizeGroup boundaries
        PtrList<dimensionedScalar> v_;

        //- Section width required for binary breakup formulation
        PtrList<PtrList<dimensionedScalar>> delta_;

        //- Explicitly treated sources
        PtrList<volScalarField> Su_;

        //- Implicitly treated sources
        PtrList<volScalarField> Sp_;

        //- Dilatation errors
        HashTable<volScalarField> dilatationErrors_;

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

        void initialiseDmdtfs();

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

        void correctDilatationError();

        void calcAlphas();

        tmp<volScalarField> calcDsm();

        void calcVelocity();

        //- Return whether the sources should be updated on this iteration
        bool updateSources();

        //- Return the interval at which the sources are updated
        inline label sourceUpdateInterval() const;

public:

    //- Runtime type information
    TypeName("populationBalanceModel");


    // Constructors

        //- Construct for a fluid
        populationBalanceModel(const phaseSystem& fluid, const word& name);

        //- Return clone
        autoPtr<populationBalanceModel> clone() const;

        //- Return a pointer to a new populationBalanceModel object created on
        //  freestore from Istream
        class iNew
        {
            const phaseSystem& fluid_;

        public:

            iNew(const phaseSystem& fluid)
            :
                fluid_(fluid)
            {}

            autoPtr<populationBalanceModel> operator()(Istream& is) const
            {
                const word name(is);

                Info << "Setting up population balance: " << name << endl;

                return autoPtr<populationBalanceModel>
                (
                    new populationBalanceModel(fluid_, name)
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

        //- Return reference to the interfacial mass transfer rates
        inline const phaseSystem::dmdtfTable& dmdtfs() const;

        //- Return reference to the mesh
        inline const fvMesh& mesh() const;

        //- Return populationBalanceCoeffs dictionary
        inline const dictionary& dict() const;

        //- Return the number of corrections
        inline label nCorr() const;

        //- Solve on final pimple iteration only
        inline Switch solveOnFinalIterOnly() const;

        //- Return continuous phase
        inline const phaseModel& continuousPhase() const;

        //- Return the velocity groups belonging to this populationBalance
        inline const HashTable<const velocityGroup*>& velocityGroupPtrs() const;

        //- Return the size groups belonging to this populationBalance
        inline const UPtrList<sizeGroup>& sizeGroups() const;

        //- Return implicit source terms
        inline const volScalarField& Sp(const label i) const;

        //- Return dilatation obtained from source terms
        inline const HashTable<volScalarField>& sourceDilatation() const;

        //- Return coalescence relevant size group pairs
        inline const List<labelPair>& coalescencePairs() const;

        //- Return binary breakup relevant size group pairs
        inline const List<labelPair>& binaryBreakupPairs() const;

        //- Return total void of phases belonging to this populationBalance
        inline const volScalarField& alphas() const;

        //- Return average velocity
        inline const volVectorField& U() const;

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

        //- Return reference to momentumTransport model of the continuous phase
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
