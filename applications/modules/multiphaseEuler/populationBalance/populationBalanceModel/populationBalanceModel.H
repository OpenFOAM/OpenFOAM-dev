/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2017-2025 OpenFOAM Foundation
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
    Foam::populationBalanceModel

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

    Verification cases for the population balance modelling functionality are
    provided in test/multiphaseEuler/populationBalance.

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
    // Define the phases
    phases          (air water);

    // Settings for the air phase. This phase is in the form of bubbles, the
    // diameters of which will be represented by a population balance model.
    // So, set the diameter model to 'populationBalance' and set the name of
    // the associated population balance to 'bubbles'. The name of the
    // population balance is arbitrary; it just needs to be unique and
    // descriptive as it will be referred to below. And if multiple phases
    // share the same population balance then it needs to be the same for all
    // those phases.
    air
    {
        type            pureIsothermalPhaseModel;

        diameterModel
        {
            type            populationBalance;
            populationBalance bubbles;
            nGroups         5;
        }

        residualAlpha   1e-6;
    }

    // Settings for the water phase. This phase is continuous, so there is no
    // diameter model.
    water
    {
        type            pureIsothermalPhaseModel;

        diameterModel   none;

        residualAlpha   1e-6;
    }

    // Set the blending to 'continuous' so that only sub-models relating to the
    // air_dispersedIn_water configuration are used
    blending
    {
        default
        {
            type            continuous;
            phase           water;
        }
    }

    // Set the surface tension
    surfaceTension
    {
        air_water
        {
            type            constant;
            sigma           0.07;
        }
    }

    // Settings for the bubbles population balance model
    bubbles
    {
        // The continuous phase surrounding the bubbles is the water phase
        continuousPhase water;

        // Give the groups uniformly spaced diameters between 1 [mm] and 5 [mm]
        sphericalDiameters
        {
            type            uniform;
            min             1e-3;
            max             5e-3;
        }

        // Define the bubbles to be spherical
        shapeModel      spherical;

        // Choose coalescence and breakup sub-models ...

        coalescenceModels
        (
            LehrMilliesMewes
            {}
        );

        breakupModels
        ();

        binaryBreakupModels
        (
            LehrMilliesMewes
            {}
        );
    }
    \endverbatim

See also
    Foam::populationBalanceSystem
    Foam::diameterModels::populationBalance
    Foam::populationBalance::coalescenceModel
    Foam::populationBalance::breakupModel
    Foam::populationBalance::daughterSizeDistributionModel
    Foam::populationBalance::binaryBreakupModel
    Foam::functionObjects::populationBalanceSizeDistribution
    Foam::functionObjects::populationBalanceSetSizeDistribution
    Foam::functionObjects::populationBalanceMoments

SourceFiles
    populationBalanceModel.C

\*---------------------------------------------------------------------------*/

#ifndef populationBalanceModel_H
#define populationBalanceModel_H

#include "populationBalance.H"
#include "phaseSystem.H"
#include "HashPtrTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class distribution;
class phaseCompressibleMomentumTransportModel;

namespace populationBalance
{
    class shapeModel;
    class coalescenceModel;
    class breakupModel;
    class binaryBreakupModel;
}

/*---------------------------------------------------------------------------*\
                   Class populationBalanceModel Declaration
\*---------------------------------------------------------------------------*/

class populationBalanceModel
:
    public regIOobject
{
public:

    // Public Type Definitions

        //- Table of interfacial mass transfer rates
        typedef
            HashPtrTable
            <
                volScalarField::Internal,
                phaseInterfaceKey,
                phaseInterfaceKey::hash
            >
            dmdtfTable;


private:

    // Private Data

        //- Reference to the phaseSystem
        const phaseSystem& fluid_;

        //- Reference to the mesh
        const fvMesh& mesh_;

        //- Name of the populationBalance
        const word name_;

        //- Continuous phase
        const phaseModel& continuousPhase_;

        //- The phase associated with each group
        UPtrList<const phaseModel> phases_;

        //- A minimal list of the phases in this population balance
        UPtrList<const phaseModel> uniquePhases_;

        //- The diameter model associated with each group
        UPtrList<const diameterModels::populationBalance> diameters_;

        //- A minimal list of the diameter models in this population balance
        UPtrList<const diameterModels::populationBalance> uniqueDiameters_;

        //- The group fractions
        PtrList<volScalarField> fs_;

        //- The representative spherical diameters of the groups
        PtrList<dimensionedScalar> dSphs_;

        //- The representative volumes of the groups
        PtrList<dimensionedScalar> vs_;

        //- Explicit coalescence and breakup source terms
        PtrList<volScalarField::Internal> Su_;

        //- Implicit coalescence and breakup source terms
        PtrList<volScalarField::Internal> Sp_;

        //- Coalescence and breakup mass transfer rates
        dmdtfTable dmdtfs_;

        //- Expansion mass transfer rates
        dmdtfTable expansionDmdtfs_;

        //- Model source mass transfer rates
        dmdtfTable modelSourceDmdtfs_;

        //- Rates of change of volume per unit volume
        PtrList<volScalarField::Internal> expansionRates_;

        //- Dilatation errors
        PtrList<volScalarField::Internal> dilatationErrors_;

        //- The shape model
        autoPtr<populationBalance::shapeModel> shapeModel_;

        //- Coalescence models
        PtrList<populationBalance::coalescenceModel> coalescenceModels_;

        //- Coalescence rate
        autoPtr<volScalarField::Internal> coalescenceRate_;

        //- Coalescence relevant group pairs
        List<labelPair> coalescencePairs_;

        //- Breakup models
        PtrList<populationBalance::breakupModel> breakupModels_;

        //- Breakup rate
        autoPtr<volScalarField::Internal> breakupRate_;

        //- Binary breakup models
        PtrList<populationBalance::binaryBreakupModel> binaryBreakupModels_;

        //- Binary breakup rate
        autoPtr<volScalarField::Internal> binaryBreakupRate_;

        //- Section width required for binary breakup formulation
        PtrList<PtrList<dimensionedScalar>> binaryBreakupDeltas_;

        //- Binary breakup relevant group pairs
        List<labelPair> binaryBreakupPairs_;

        //- Total void fraction
        autoPtr<volScalarField> alphas_;

        //- Mean Sauter diameter
        autoPtr<volScalarField> dsm_;

        //- Average velocity
        autoPtr<volVectorField> U_;

        //- Counter for interval between source term updates
        label sourceUpdateCounter_;


    // Private Member Functions

        //- Get the coefficient dictionary
        const dictionary& coeffDict() const;

        //- Do pre-calculations relating to coalescence and breakup
        void precomputeCoalescenceAndBreakup();

        //- Calculate particle creation by coalescence
        void birthByCoalescence(const label j, const label k);

        //- Calculate particle destruction by coalescence
        void deathByCoalescence(const label i, const label j);

        //- Calculate particle creation by breakup
        void birthByBreakup(const label k, const label model);

        //- Calculate particle destruction by breakup
        void deathByBreakup(const label i);

        //- Calculate particle creation by binary-breakup
        void birthByBinaryBreakup(const label i, const label j);

        //- Calculate particle destruction by binary-breakup
        void deathByBinaryBreakup(const label j, const label i);

        //- Do calculations relating to coalescence and breakup
        void computeCoalescenceAndBreakup();

        //- Do pre-calculations relating to expansion
        void precomputeExpansion();

        //- Return the expansion sources for a group fraction, or (if flds is
        //  valid) a secondary property
        Pair<tmp<volScalarField::Internal>> expansionSus
        (
            const label i,
            const UPtrList<volScalarField>& flds = UPtrList<volScalarField>()
        ) const;

        //- Do calculations relating to expansion
        void computeExpansion();

        //- Do pre-calculations relating to model sources
        void precomputeModelSources();

        //- Return the model sources for a group fraction, or (if flds is
        //  valid) a secondary property
        Pair<tmp<volScalarField::Internal>> modelSourceRhoSus
        (
            const label i,
            const UPtrList<volScalarField>& flds = UPtrList<volScalarField>()
        ) const;

        //- Do calculations relating to model sources
        void computeModelSources();

        //- Update the dilatation errors
        void computeDilatationErrors();

        //- Return whether the sources should be updated on this iteration
        bool updateSources();

        //- Return the interval at which the sources are updated
        inline label sourceUpdateInterval() const;

        //- Return the coefficients of the lower half of the number allocation
        //  coefficient function. This spans the range from group i-1 to group
        //  i. The first coefficient is a constant, and the second is
        //  multiplied by v.
        Pair<dimensionedScalar> etaCoeffs0(const label i) const;

        //- Return the coefficients of the lower half of the number allocation
        //  coefficient function. This spans the range from group i to group
        //  i+1. The first coefficient is a constant, and the second is
        //  multiplied by v.
        Pair<dimensionedScalar> etaCoeffs1(const label i) const;

        //- Return the coefficients of the lower half of the volume allocation
        //  coefficient function. This spans the range from group i-1 to group
        //  i. The first coefficient is a constant, and the second is
        //  multiplied by 1/v.
        Pair<dimensionedScalar> etaVCoeffs0(const label i) const;

        //- Return the coefficients of the lower half of the volume allocation
        //  coefficient function. This spans the range from group i to group
        //  i+1. The first coefficient is a constant, and the second is
        //  multiplied by 1/v.
        Pair<dimensionedScalar> etaVCoeffs1(const label i) const;


public:

    //- Runtime type information
    TypeName("populationBalanceModel");


    // Static Member Functions

        //- Return IO for a group-associated field
        static IOobject groupFieldIo
        (
            const word& name,
            const label i,
            const phaseModel& phase,
            const IOobject::readOption r = IOobject::NO_READ,
            const bool registerObject = true
        );

        //- Read and return a group-associated field
        static tmp<volScalarField> groupField
        (
            const word& name,
            const label i,
            const phaseModel& phase
        );


    // Constructors

        //- Construct for a fluid and with a given name
        populationBalanceModel(const phaseSystem& fluid, const word& name);


    //- Destructor
    virtual ~populationBalanceModel();


    // Member Functions

        //- Return reference to the phaseSystem
        inline const phaseSystem& fluid() const;

        //- Return reference to the mesh
        inline const fvMesh& mesh() const;

        //- Return continuous phase
        inline const phaseModel& continuousPhase() const;

        //- Return the number of groups in the population balance
        inline label nGroups() const;

        //- Access the list of phases associated with each group
        inline const UPtrList<const phaseModel>& phases() const;

        //- Access the minimal list of the phases in this population balance
        inline const UPtrList<const phaseModel>& uniquePhases() const;

        //- Access the list of diameter models associated with each group
        inline const UPtrList<const diameterModels::populationBalance>&
            diameters() const;

        //- Access the minimal list of the diameter models in this population
        //  balance
        inline const UPtrList<const diameterModels::populationBalance>&
            uniqueDiameters() const;

        //- Access the group fractions
        inline const PtrList<volScalarField>& fs() const;

        //- Modify the group fractions
        inline PtrList<volScalarField>& fs();

        //- Access a group fraction
        inline const volScalarField& f(const label i) const;

        //- Modify a group fraction
        inline volScalarField& f(const label i);

        //- Access the representative spherical diameters of the groups
        inline const PtrList<dimensionedScalar>& dSphs() const;

        //- Access the representative spherical diameters of a group
        inline const dimensionedScalar& dSph(const label i) const;

        //- Access the representative volumes diameters of the groups
        inline const PtrList<dimensionedScalar>& vs() const;

        //- Access the representative volumes diameters of a group
        inline const dimensionedScalar& v(const label i) const;

        //- Access the shape model
        const populationBalance::shapeModel& shape() const;

        //- Return the representative surface area of the group
        tmp<volScalarField> a(const label i) const;

        //- Return the representative diameter of the group
        tmp<volScalarField> d(const label i) const;

        //- Return reference to the coalescence and breakup interfacial mass
        //  transfer rates
        inline const dmdtfTable& dmdtfs() const;

        //- Return reference to the expansion interfacial mass transfer rates
        inline const dmdtfTable& expansionDmdtfs() const;

        //- Return reference to the model source interfacial mass transfer rates
        inline const dmdtfTable& modelSourceDmdtfs() const;

        //- Return solution settings dictionary for this populationBalance
        inline const dictionary& solverDict() const;

        //- Solve on final pimple iteration only
        inline bool solveOnFinalIterOnly() const;

        //- Return coalescence relevant group pairs
        inline const List<labelPair>& coalescencePairs() const;

        //- Return binary breakup relevant group pairs
        inline const List<labelPair>& binaryBreakupPairs() const;

        //- Return total void of phases belonging to this populationBalance
        inline const volScalarField& alphas() const;

        //- Return average velocity
        inline tmp<volVectorField> U() const;

        //- Return the number allocation coefficient for a single volume
        dimensionedScalar eta(const label i, const dimensionedScalar& v) const;

        //- Return the number allocation coefficient for a field of volumes
        tmp<volScalarField::Internal> eta
        (
            const label i,
            const volScalarField::Internal& v
        ) const;

        //- Return the volume allocation coefficient for a single volume
        dimensionedScalar etaV(const label i, const dimensionedScalar& v) const;

        //- Return the volume allocation coefficient for a field of volumes
        tmp<volScalarField::Internal> etaV
        (
            const label i,
            const volScalarField::Internal& v
        ) const;

        //- Return the volume allocation coefficient for a given distribution
        //  of diameters over a range of groups. The diameter distribution
        //  should be volumetrically sampled (i.e., sampleQ should equal 3).
        dimensionedScalar etaV(const labelPair is, const distribution& d) const;

        //- Return the volume allocation coefficient for a given distribution
        //  of diameters. The diameter distribution should be volumetrically
        //  sampled (i.e., sampleQ should equal 3). The result is normalised
        //  with respect to the volume allocation coefficient for the phase.
        dimensionedScalar etaV(const label i, const distribution& d) const;

        //- Return the surface tension coefficient between a given dispersed
        //  and the continuous phase
        const tmp<volScalarField> sigmaWithContinuousPhase
        (
            const phaseModel& dispersedPhase
        ) const;

        //- Return the surface tension coefficient between a given group
        //  and the continuous phase
        const tmp<volScalarField> sigmaWithContinuousPhase(const label i) const;

        //- Return reference to momentumTransport model of the continuous phase
        const phaseCompressibleMomentumTransportModel&
            continuousTurbulence() const;

        //- Return the implicit coalescence and breakup source term
        tmp<volScalarField::Internal> Sp(const label i) const;

        //- Return the explicit expansion source term
        tmp<volScalarField::Internal> expansionSu
        (
            const label i,
            const UPtrList<volScalarField>& flds = UPtrList<volScalarField>()
        ) const;

        //- Return the implicit expansion source term
        tmp<volScalarField::Internal> expansionSp(const label i) const;

        //- Return the explicit model source source term
        tmp<volScalarField::Internal> modelSourceSu
        (
            const label i,
            const UPtrList<volScalarField>& flds = UPtrList<volScalarField>()
        ) const;

        //- Solve the population balance equation
        void solve();

        //- Correct derived quantities
        void correct();

        //- Dummy write for regIOobject
        bool writeData(Ostream&) const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "populationBalanceModelI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
