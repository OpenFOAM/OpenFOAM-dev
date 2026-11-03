/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2026 OpenFOAM Foundation
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
    Foam::turbulentSyntheticEddyInletFvPatchVectorField

Description
    This boundary condition generates a synthetic-turbulence velocity inlet
    using a divergence-free synthetic eddy method. Fluctuations are
    constructed from a population of virtual eddies convected across the
    patch, reproducing a target Reynolds stress tensor and integral length
    scale field.

    By default the target profiles are obtained from bundled DNS channel
    flow data (Re_tau = 180, 395 and 590), mapped onto the patch in wall
    units. Alternatively, user-supplied fields of the Reynolds stress
    (volSymmTensorField), integral length scale (volScalarField) and mean
    velocity (volVectorField) may be referenced by name.

    Wall and cyclic neighbour patches are detected automatically from the
    mesh, eddy spawn positions are drawn by area-weighted face sampling,
    and swirl is generalised through a user-overridable swirl axis and
    centre.

    Two optional self-corrections are provided: an initialisation-time
    Monte-Carlo correction of the inlet Reynolds stresses, and a runtime
    correction which compensates for decay between the inlet and a
    downstream sampling plane.

    The inflow magnitude is specified either as a mean velocity or as a
    mass flow rate, both of which may be time-varying Function1 entries.

    Reference:
    \verbatim
        Fyhn, E. H., & Fyhn, H. (2026).
        turbulentSyntheticEddyInlet: A self-correcting, divergence-free
        synthetic-eddy boundary condition for OpenFOAM.
        SoftwareX, 35, 102907.
    \endverbatim

Usage
    \table
        Property     | Description                  | Required | Default value
        Umean        | mean inlet velocity          | if massFlow absent |
        massFlow     | mass flow rate (alternative to Umean) | no |
        rhoInf       | reference density for massFlow | no     | 1
        Utau         | friction velocity for wall-unit scaling | no | estimated
        wallPatches  | wall patches for wall distance | no     | detected
        eddyDensity  | eddy number density factor   | no       | 5
        eddyConvectionFactor | eddy convection speed factor | no | 0.5
        nCellPerEddy | minimum cells per eddy       | no       | 3
        Lscale       | length scale multiplier      | no       | 1
        Rscale       | Reynolds stress multiplier   | no       | 1
        swirlNumber  | inlet swirl number           | no       | 0
        swirlAxis    | swirl axis                   | no       | detected
        swirlCenter  | swirl centre                 | no       | detected
        customRName  | name of user Reynolds stress field | no | DNS profile
        customLName  | name of user length scale field | no    | DNS profile
        customUmeanName | name of user mean velocity field | no | DNS profile
        inletRcorrection | Monte-Carlo inlet stress correction | no | false
        inletRcorrectionOverlaps | Monte-Carlo eddy-overlap density | no | 500
        correctionDistance | downstream correction distance; 0 = off | no | 0
        downstreamStartTime | start time of downstream correction | no | 0
        downstreamMaxPasses | maximum correction passes | no    | 10
        downstreamWindowEddies | averaging window, in eddy passes | no | 10
        downstreamMaxAmp | maximum correction amplification | no | 10
        downstreamPatience | non-improving passes before stopping | no | 3
        downstreamTol | relative convergence tolerance | no     | 0.05
        seed         | random number seed; -1 = randomise | no  | 536
    \endtable

    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            turbulentSyntheticEddyInlet;
        massFlow        2.53e-4;
    }
    \endverbatim

See also
    Foam::fixedValueFvPatchField
    Foam::turbulentInletFvPatchField

SourceFiles
    turbulentSyntheticEddyInletFvPatchVectorField.C
    syntheticEddy.C

\*---------------------------------------------------------------------------*/

#ifndef turbulentSyntheticEddyInletFvPatchVectorField_H
#define turbulentSyntheticEddyInletFvPatchVectorField_H

#include "fixedValueFvPatchFields.H"
#include "syntheticEddy.H"
#include "randomGenerator.H"
#include "Function1.H"
#include "dnsCatalog.H"


namespace Foam
{

class turbulentSyntheticEddyInletFvPatchVectorField
:
    public fixedValueFvPatchVectorField
{
    // Friction velocity and y -> y+ scaling
    scalar Utau_;
    bool UtauIsInDict_;
    scalar yPlusScale_;
    scalar avNu_;

    // Mean velocity at the inlet (unit-mean profile)
    vectorField U_;

    // Inflow specification
    bool massFlowMode_;
    scalar eddyConvectionFactor_;

    autoPtr<Function1<scalar>> inFlow_;
    scalar Umean_;
    scalar rhoInf_;

    // Patch geometry (auto-detected)
    vector centerCoord_;
    vector patchNormal_;
    scalar delta_;

    // Walls used for wall distance and local-frame rotation
    wordList wallPatchNames_;
    bool wallPatchesIsInDict_;
    scalarField wallDist_;
    vectorField wallNormal_;
    scalar maxVertexWallDist_;

    // Swirl
    scalar swirlNumber_;
    vector swirlAxis_;
    vector swirlCenter_;
    bool swirlAxisIsInDict_;
    bool swirlCenterIsInDict_;

    // Per-face precomputed eddy data (local to processor)
    scalarField faceSigmaX_;
    tensorField faceRpg_;
    scalarField faceGamma_;
    vectorField faceAlphaCoeff_;
    boolList faceActive_;
    symmTensorField faceRlocal_;
    scalarField faceLlocal_;

    // User-supplied R/L bypass the Utau^2 scaling at spawn time
    bool useUserR_;
    bool useUserL_;

    // Optional names of registered fields (volSymmTensorField/volScalarField/
    // volVectorField) overriding the DNS R/L/Umean profile on this patch.
    // Empty -> use the bundled DNS profile.
    word customRName_;
    word customLName_;
    word customUmeanName_;

    // Init-time MC pre-distortion of the target R so the BC's spatial-averaging
    // output matches the prescribed target R.
    bool inletRcorrection_;

    // MC eddy-overlap density: N = clamp(mcOverlaps_*V_box/V_eddy, 1000, 5e6).
    scalar mcOverlaps_;

    // Globally-gathered face data
    pointField gCf_;
    scalarField gMagSf_;
    scalarField gSigmaX_;
    tensorField gRpg_;
    scalarField gGamma_;
    vectorField gAlphaCoeff_;
    boolList gActive_;
    scalarField gCDF_;

    // Spanwise cyclic ghost translations (periodic eddy images)
    List<vector> cyclicTranslations_;

    // Eddy state
    List<syntheticEddy> eddies_;
    label nEddy_;
    const scalar eddyDensity_;
    const scalar Lscale_;
    const scalar Rscale_;
    const label nCellPerEddy_;
    scalar v0_;
    scalar maxSigmaX_;

    // DNS dataset chosen at init time (closest Re_tau)
    const DNSDataset* selectedDataset_;

    // RNG
    // seedInput_ is the raw dict value (seed = -1 -> re-randomised on every init)
    const label seedInput_;
    const label seed_;
    randomGenerator randGen_;

    // First-call deferred init
    bool initialized_;

    // Helpers
    void updateUmean();
    //- Friction velocity from the bulk velocity via a skin-friction relation.
    scalar frictionUtau(scalar Ub, scalar nu) const;
    tmp<scalarField> rhoOnPatch() const;
    void initializeGeometry();
    void computeWallDistance();
    void initializeProfiles();
    void buildGhostTransforms();
    void gatherFaceData();

    // Per-face buildEddyTemplate loop from the current faceRlocal_/faceLlocal_.
    void rebuildTemplates();
    // Init-time Monte-Carlo iterative pre-distortion of eddy intensities.
    void applyInletRcorrectionMC();

    // Shared sample set for the R-correction
    struct RcSampleSet
    {
        bool valid = false;
        // Global active faces.
        label nF = 0;
        scalar A_total = 0;
        List<vector> CfF, wnF;
        scalarList wdF, aF, LF;
        // Per-face derived.
        scalar h = 0, cellDxAvg = 0, maxSigma = 0;
        scalar wnExcludeCoeff = 0.15;   // near-wall exclusion band coeff
        Field<symmTensor> R_ref_face;
        Field<tensor> Tframe_face;
        scalarField sigmaX_face;
        // Samples (sparse subset).
        label nS = 0;
        List<vector> Cf_s;
        scalarList a_s;
        Field<symmTensor> R_ref_s;
        Field<tensor> Tframe_s;
        scalar sigmaKernel = 0;
        // Kernel maps.
        List<labelList> nbrs_ss, nbrs_fs;
        List<scalarField> W_ss, W_fs;
        scalarField wnorm_ss, wnorm_fs;
        // Write-back of a face-space field to local faceRlocal_ entries.
        labelList origIdxLocal;
        label off = 0;
    };
    // Gather + Poisson-disk selection + kernel-map build.
    void buildSampleSet(RcSampleSet& ss) const;
    // Smooth a sample-space symmTensor field onto every global face.
    void scatterSampleToFace
    (
        const RcSampleSet& ss,
        const Field<symmTensor>& Rs,
        Field<symmTensor>& Rface
    ) const;
    // Area-weighted, component-normalised RMS residual.
    scalar rmsResidual(const RcSampleSet& ss, const Field<symmTensor>& diff) const;

    //- Eddy u' at xp including all cyclic ghost images
    vector eddyImageContribution
    (
        const syntheticEddy& ed, const point& pos, scalar reach, const point& xp
    ) const;

    //- Spawn nEddies and accumulate their mean-removed R into Rsum (sample space).
    void spawnEddyBatch
    (
        const RcSampleSet& ss, const scalarList& cumA, randomGenerator& rng,
        label nEddies, const Field<symmTensor>& R_target_face,
        symmTensorField& Rsum
    ) const;

    //- Debug: write a face-local-frame R as a global-frame volSymmTensorField.
    void writeRcDebugField
    (
        const word& name, const Field<symmTensor>& faceLocalDict,
        const RcSampleSet& ss, scalar physScale
    ) const;

    // ---- Downstream-measurement R-correction ----
    struct downstreamSettings
    {
        scalar distance = 0;        // correctionDistance [m], 0 disables
        scalar windowEddies = 10;   // window = eddies*tau_eddy (tau_eddy = maxSigmaX/U_b)
        label  maxPasses = 10;      // cap on closed-loop passes
        scalar startTime = 0;       // gate before pass 1 (case develops)
        scalar maxAmp = 10;         // cap the global scale so no normal stress exceeds maxAmp*reference
        label  patience = 3;        // stop after this many non-improving passes
        scalar tol = 0.05;          // converged when TKE error stays below this
    };
    downstreamSettings ds_;

    // Runtime controller state (transient, reset on any clone/map).
    struct downstreamState
    {
        enum class phase { warmup, accumulate, done };
        bool setupDone = false;         // probe ownership built?
        phase ph = phase::warmup;
        label pass = 0;
        scalar phaseStartTime = 0;      // sim time the current phase began
        scalar sampleCount = 0;         // timesteps accumulated this window
        scalar windowEddies = 0;        // current window size (grows on stall)
        scalar maxSigma = 0;            // global max sigmaX (tau_eddy = sigmaX/U_b)

        // Probe arrays: one probe per active inlet face
        pointField probeX;             // Cf + distance*patchNormal
        scalarField probeWd;           // wallDist of the source face (diagnostics)
        scalarField probeArea;         // face area (mean weighting)
        vectorField probeWn;           // wallNormal (local-frame rotation)
        labelList ownerProc;           // min-rank owner per probe
        labelList localCell;           // local celli (this rank) or -1
        labelList myProbes;            // probe indices this rank samples
        vectorField sumU;              // sum U accumulator (per probe)
        symmTensorField sumUU;         // sum U(x)U accumulator (per probe)

        // Reference profile
        Field<symmTensor> probeRef;    // per-probe reference R (local frame)
        Field<symmTensor> probeTarget; // per-probe applied inlet target (diag)
        symmTensor refMean = Zero;     // area-weighted mean reference R
        symmTensor appliedMean = Zero; // area-weighted mean of the applied target

        symmTensorField bestTarget;    // faceRlocal_ at the lowest TKE error seen
        scalar bestTkeErr = GREAT;
        label noImprove = 0;           // consecutive non-improving passes
        label belowTol = 0;            // consecutive passes with TKE error below tol
    };
    downstreamState dsState_;

    // Reference target per local face (bundled DNS or user-supplied R)
    symmTensorField faceRlocalRef_;

    // One-time probe placement + parallel ownership assignment.
    void setupDownstreamProbes();
    // Per-timestep state machine, called from updateCoeffs after updateUmean.
    void stepDownstreamController();
    // Close an averaging window: reduce sums, form measured R, apply feedback.
    void finalizeDownstreamWindow();
    // Rebuild templates and refresh the gathered arrays after faceRlocal_ changes.
    void applyDownstreamTarget();
    // Write the converged inlet R to the startTime folder as "R_converged"
    // for restart via customRName.
    void writeConvergedR() const;

    syntheticEddy oneEddy(scalar depth);


public:

    // Per-face synthetic-eddy template: time-constant shape and amplitude
    struct EddyTemplate
    {
        bool active = false;      // false if R violates admissibility
        tensor Rpg = tensor::I;   // principal-to-global rotation
        vector alphaCoeff = Zero; // unsigned alpha magnitudes
        scalar gamma = 1;         // aspect ratio sigma_x / sigma_y
        scalar sigmaX = 0;        // major (streamwise) length scale
        scalar c1 = -1;           // normalisation factor
    };

    // Build an eddy template from a local-frame R and L_h. 
    // The local frame is wall-aligned: u = patchNormal, v = wallNormal, w = u x v.
    static EddyTemplate buildEddyTemplate
    (
        const symmTensor& Rlocal,
        scalar Lh,
        const vector& patchNormal,
        const vector& wallNormal,
        scalar Lscale,
        scalar halfWidth,
        scalar cellDx,
        label nCellPerEddy,
        bool useUserR = false,
        bool useUserL = false
    );

private:

    //- Resolve the RNG seed. Returns the request unchanged, except seed = -1
    //  draws a random seed on the master and broadcasts it to all ranks.
    static label resolveSeed(label requested);


public:

    TypeName("turbulentSyntheticEddyInlet");


    turbulentSyntheticEddyInletFvPatchVectorField
    (
        const fvPatch&,
        const DimensionedField<vector, fvMesh>&,
        const dictionary&
    );

    turbulentSyntheticEddyInletFvPatchVectorField
    (
        const turbulentSyntheticEddyInletFvPatchVectorField&,
        const fvPatch&,
        const DimensionedField<vector, fvMesh>&,
        const fieldMapper&
    );

    turbulentSyntheticEddyInletFvPatchVectorField
    (
        const turbulentSyntheticEddyInletFvPatchVectorField&
    ) = delete;

    turbulentSyntheticEddyInletFvPatchVectorField
    (
        const turbulentSyntheticEddyInletFvPatchVectorField&,
        const DimensionedField<vector, fvMesh>&
    );

    virtual tmp<fvPatchVectorField> clone
    (
        const DimensionedField<vector, fvMesh>& iF
    ) const
    {
        return tmp<fvPatchVectorField>
        (
            new turbulentSyntheticEddyInletFvPatchVectorField(*this, iF)
        );
    }


private:

    void initializeEddies();
    void convectEddies();
    vector uPrimeEddy(const vector& pos) const;


public:

    virtual void updateCoeffs();
    virtual void write(Ostream&) const;
};

} // End namespace Foam

#endif

// ************************************************************************* //
