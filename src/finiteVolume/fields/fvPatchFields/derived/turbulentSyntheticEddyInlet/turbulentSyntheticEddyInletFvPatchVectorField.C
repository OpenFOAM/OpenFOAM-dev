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

\*---------------------------------------------------------------------------*/

#include "turbulentSyntheticEddyInletFvPatchVectorField.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "fvcMeshPhi.H"
#include "addToRunTimeSelectionTable.H"
#include "mathematicalConstants.H"
#include "wallPolyPatch.H"
#include "cyclicPolyPatch.H"
#include "Switch.H"
#include "Pstream.H"
#include "scalarMatrices.H"
#include "meshSearch.H"
#include "clock.H"
#include "clockTime.H"
#include "ListOps.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::turbulentSyntheticEddyInletFvPatchVectorField::
turbulentSyntheticEddyInletFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, fvMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchVectorField(p, iF, dict, false),
    Utau_(dict.lookupOrDefault<scalar>("Utau", 0.0)),
    UtauIsInDict_(dict.found("Utau", false, false)),
    yPlusScale_(0),
    avNu_(0),
    U_(size(), Zero),
    massFlowMode_(dict.found("massFlow")),
    eddyConvectionFactor_
    (
        dict.lookupOrDefault<scalar>("eddyConvectionFactor", 0.5)
    ),
    inFlow_
    (
        Function1<scalar>::New
        (
            massFlowMode_ ? "massFlow" : "Umean",
            this->db().time().userUnits(),
            iF.dimensions(),
            dict
        )
    ),
    Umean_(0),
    rhoInf_(dict.lookupOrDefault<scalar>("rhoInf", 1)),
    centerCoord_(Zero),
    patchNormal_(Zero),
    delta_(1),
    wallPatchNames_(dict.lookupOrDefault<wordList>("wallPatches", wordList())),
    wallPatchesIsInDict_(dict.found("wallPatches")),
    maxVertexWallDist_(0),
    swirlNumber_(dict.lookupOrDefault<scalar>("swirlNumber", 0)),
    swirlAxis_(dict.lookupOrDefault<vector>("swirlAxis", Zero)),
    swirlCenter_(dict.lookupOrDefault<vector>("swirlCenter", Zero)),
    swirlAxisIsInDict_(dict.found("swirlAxis")),
    swirlCenterIsInDict_(dict.found("swirlCenter")),
    useUserR_(false),
    useUserL_(false),
    customRName_(dict.lookupOrDefault<word>("customRName", word())),
    customLName_(dict.lookupOrDefault<word>("customLName", word())),
    customUmeanName_(dict.lookupOrDefault<word>("customUmeanName", word())),
    inletRcorrection_(dict.lookupOrDefault<Switch>("inletRcorrection", false)),
    mcOverlaps_(dict.lookupOrDefault<scalar>("inletRcorrectionOverlaps", 500)),
    nEddy_(0),
    eddyDensity_(dict.lookupOrDefault<scalar>("eddyDensity", 5)),
    Lscale_(dict.lookupOrDefault<scalar>("Lscale", 1)),
    Rscale_(dict.lookupOrDefault<scalar>("Rscale", 1)),
    nCellPerEddy_(dict.lookupOrDefault<label>("nCellPerEddy", 3)),
    v0_(0),
    maxSigmaX_(0),
    selectedDataset_(nullptr),
    seedInput_(dict.lookupOrDefault<label>("seed", 536)),
    seed_(resolveSeed(seedInput_)),
    randGen_(seed_),
    initialized_(false)
{
    // Downstream-correction settings
    ds_.distance      = dict.lookupOrDefault<scalar>("correctionDistance", 0);
    ds_.windowEddies  = dict.lookupOrDefault<scalar>("downstreamWindowEddies", 10);
    ds_.maxPasses     = dict.lookupOrDefault<label>("downstreamMaxPasses", 10);
    ds_.startTime     = dict.lookupOrDefault<scalar>("downstreamStartTime", 0);
    ds_.maxAmp        = dict.lookupOrDefault<scalar>("downstreamMaxAmp", 10);
    ds_.patience      = dict.lookupOrDefault<label>("downstreamPatience", 3);
    ds_.tol           = dict.lookupOrDefault<scalar>("downstreamTol", 0.05);

    vectorField::operator=(U_);
}


Foam::turbulentSyntheticEddyInletFvPatchVectorField::
turbulentSyntheticEddyInletFvPatchVectorField
(
    const turbulentSyntheticEddyInletFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, fvMesh>& iF,
    const fieldMapper& mapper
)
:
    fixedValueFvPatchVectorField(ptf, p, iF, mapper, false),
    Utau_(ptf.Utau_),
    UtauIsInDict_(ptf.UtauIsInDict_),
    yPlusScale_(ptf.yPlusScale_),
    avNu_(ptf.avNu_),
    U_(ptf.U_),
    massFlowMode_(ptf.massFlowMode_),
    eddyConvectionFactor_(ptf.eddyConvectionFactor_),
    inFlow_(ptf.inFlow_, false),
    Umean_(ptf.Umean_),
    rhoInf_(ptf.rhoInf_),
    centerCoord_(ptf.centerCoord_),
    patchNormal_(ptf.patchNormal_),
    delta_(ptf.delta_),
    wallPatchNames_(ptf.wallPatchNames_),
    wallPatchesIsInDict_(ptf.wallPatchesIsInDict_),
    maxVertexWallDist_(ptf.maxVertexWallDist_),
    swirlNumber_(ptf.swirlNumber_),
    swirlAxis_(ptf.swirlAxis_),
    swirlCenter_(ptf.swirlCenter_),
    swirlAxisIsInDict_(ptf.swirlAxisIsInDict_),
    swirlCenterIsInDict_(ptf.swirlCenterIsInDict_),
    useUserR_(ptf.useUserR_),
    useUserL_(ptf.useUserL_),
    customRName_(ptf.customRName_),
    customLName_(ptf.customLName_),
    customUmeanName_(ptf.customUmeanName_),
    inletRcorrection_(ptf.inletRcorrection_),
    mcOverlaps_(ptf.mcOverlaps_),
    eddies_(ptf.eddies_),
    nEddy_(ptf.nEddy_),
    eddyDensity_(ptf.eddyDensity_),
    Lscale_(ptf.Lscale_),
    Rscale_(ptf.Rscale_),
    nCellPerEddy_(ptf.nCellPerEddy_),
    v0_(ptf.v0_),
    maxSigmaX_(ptf.maxSigmaX_),
    selectedDataset_(ptf.selectedDataset_),
    seedInput_(ptf.seedInput_),
    seed_(ptf.seed_),
    randGen_(ptf.randGen_),
    initialized_(ptf.initialized_),
    ds_(ptf.ds_)
{}


Foam::turbulentSyntheticEddyInletFvPatchVectorField::
turbulentSyntheticEddyInletFvPatchVectorField
(
    const turbulentSyntheticEddyInletFvPatchVectorField& ptf,
    const DimensionedField<vector, fvMesh>& iF
)
:
    fixedValueFvPatchVectorField(ptf, iF),
    Utau_(ptf.Utau_),
    UtauIsInDict_(ptf.UtauIsInDict_),
    yPlusScale_(ptf.yPlusScale_),
    avNu_(ptf.avNu_),
    U_(ptf.U_),
    massFlowMode_(ptf.massFlowMode_),
    eddyConvectionFactor_(ptf.eddyConvectionFactor_),
    inFlow_(ptf.inFlow_, false),
    Umean_(ptf.Umean_),
    rhoInf_(ptf.rhoInf_),
    centerCoord_(ptf.centerCoord_),
    patchNormal_(ptf.patchNormal_),
    delta_(ptf.delta_),
    wallPatchNames_(ptf.wallPatchNames_),
    wallPatchesIsInDict_(ptf.wallPatchesIsInDict_),
    maxVertexWallDist_(ptf.maxVertexWallDist_),
    swirlNumber_(ptf.swirlNumber_),
    swirlAxis_(ptf.swirlAxis_),
    swirlCenter_(ptf.swirlCenter_),
    swirlAxisIsInDict_(ptf.swirlAxisIsInDict_),
    swirlCenterIsInDict_(ptf.swirlCenterIsInDict_),
    useUserR_(ptf.useUserR_),
    useUserL_(ptf.useUserL_),
    customRName_(ptf.customRName_),
    customLName_(ptf.customLName_),
    customUmeanName_(ptf.customUmeanName_),
    inletRcorrection_(ptf.inletRcorrection_),
    mcOverlaps_(ptf.mcOverlaps_),
    eddies_(ptf.eddies_),
    nEddy_(ptf.nEddy_),
    eddyDensity_(ptf.eddyDensity_),
    Lscale_(ptf.Lscale_),
    Rscale_(ptf.Rscale_),
    nCellPerEddy_(ptf.nCellPerEddy_),
    v0_(ptf.v0_),
    maxSigmaX_(ptf.maxSigmaX_),
    selectedDataset_(ptf.selectedDataset_),
    seedInput_(ptf.seedInput_),
    seed_(ptf.seed_),
    randGen_(ptf.randGen_),
    initialized_(ptf.initialized_),
    ds_(ptf.ds_)
{}


// * * * * * * * * * * * * * * * Helpers * * * * * * * * * * * * * * * * * * //

namespace
{
    using namespace Foam;
    using constant::mathematical::pi;

    // Wall-aligned local frame T (rows = u, v, w).
    inline tensor localFrameT
    (
        const vector& patchNormal, const vector& wallNormal
    )
    {
        vector u = patchNormal;
        vector v = wallNormal;
        v -= (v & u)*u;
        if (mag(v) < SMALL)
        {
            v = vector(0, 1, 0) ^ u;
            if (mag(v) < 1e-5) v = vector(0, 0, 1) ^ u;
        }
        v /= mag(v);
        const vector w = u ^ v;
        return tensor(u, v, w);
    }

    inline symmTensor globalToLocal
    (
        const symmTensor& Rglobal, const tensor& T
    )
    {
        const tensor Rt = T & tensor(Rglobal) & T.T();
        return symmTensor
        (
            Rt.xx(), Rt.xy(), Rt.xz(),
                     Rt.yy(), Rt.yz(),
                              Rt.zz()
        );
    }

    // Rotate a wall-aligned local-frame tensor back to the global mesh frame.
    inline symmTensor localToGlobal
    (
        const symmTensor& Rlocal, const tensor& T
    )
    {
        const tensor Rt = T.T() & tensor(Rlocal) & T;
        return symmTensor
        (
            Rt.xx(), Rt.xy(), Rt.xz(),
                     Rt.yy(), Rt.yz(),
                              Rt.zz()
        );
    }

    // Jacobi eigendecomposition for a 3x3 symmetric matrix (lambda
    // ascending, V columns = eigenvectors). Alternative to OpenFOAM's 
    // analytic routine, which trips on near-degenerate R.
    inline void jacobiEigen
    (
        const symmTensor& R_in, vector& lambda, tensor& V
    )
    {
        scalar A[3][3] =
        {
            {R_in.xx(), R_in.xy(), R_in.xz()},
            {R_in.xy(), R_in.yy(), R_in.yz()},
            {R_in.xz(), R_in.yz(), R_in.zz()}
        };
        scalar Q[3][3] = {{1,0,0},{0,1,0},{0,0,1}};
        for (int sweep = 0; sweep < 50; ++sweep)
        {
            const scalar off = mag(A[0][1]) + mag(A[0][2]) + mag(A[1][2]);
            if (off < 1e-14) break;
            for (int p = 0; p < 2; ++p)
                for (int q = p + 1; q < 3; ++q)
                {
                    const scalar apq = A[p][q];
                    if (mag(apq) < 1e-15) continue;
                    const scalar app = A[p][p];
                    const scalar aqq = A[q][q];
                    const scalar theta = (aqq - app)/(2*apq);
                    scalar t;
                    if (mag(theta) > 1e10) t = 0.5/theta;
                    else if (theta >= 0)
                        t = 1.0/(theta + Foam::sqrt(1.0 + theta*theta));
                    else
                        t = -1.0/(-theta + Foam::sqrt(1.0 + theta*theta));
                    const scalar c = 1.0/Foam::sqrt(1.0 + t*t);
                    const scalar s = t*c;
                    A[p][p] = app - t*apq;
                    A[q][q] = aqq + t*apq;
                    A[p][q] = A[q][p] = 0;
                    for (int r = 0; r < 3; ++r)
                    {
                        if (r == p || r == q) continue;
                        const scalar arp = A[r][p];
                        const scalar arq = A[r][q];
                        A[r][p] = A[p][r] = c*arp - s*arq;
                        A[r][q] = A[q][r] = s*arp + c*arq;
                    }
                    for (int r = 0; r < 3; ++r)
                    {
                        const scalar qrp = Q[r][p];
                        const scalar qrq = Q[r][q];
                        Q[r][p] = c*qrp - s*qrq;
                        Q[r][q] = s*qrp + c*qrq;
                    }
                }
        }
        int idx[3] = {0, 1, 2};
        const scalar lam[3] = {A[0][0], A[1][1], A[2][2]};
        for (int i = 0; i < 3; ++i)
            for (int j = i + 1; j < 3; ++j)
                if (lam[idx[i]] > lam[idx[j]])
                    Swap(idx[i], idx[j]);
        lambda = vector(lam[idx[0]], lam[idx[1]], lam[idx[2]]);
        V = tensor
        (
            Q[0][idx[0]], Q[0][idx[1]], Q[0][idx[2]],
            Q[1][idx[0]], Q[1][idx[1]], Q[1][idx[2]],
            Q[2][idx[0]], Q[2][idx[1]], Q[2][idx[2]]
        );
    }

    // Floor eigenvalues at a small fraction of the largest to keep R strictly
    // inside the PD cone, then reassemble.
    inline symmTensor projectPDjac(const symmTensor& R)
    {
        for (direction d = 0; d < 6; ++d)
        {
            if (!(mag(R[d]) < GREAT)) return symmTensor::zero;
        }
        vector lambda(Zero);
        tensor V(tensor::I);
        jacobiEigen(R, lambda, V);
        const scalar epsPD = 1e-4;               // relative PD floor
        const scalar fl = epsPD*max(lambda[2], scalar(0));  // lambda[2] = largest
        const vector lc
        (
            max(lambda[0], fl),
            max(lambda[1], fl),
            max(lambda[2], fl)
        );
        const tensor D(lc[0],0,0, 0,lc[1],0, 0,0,lc[2]);
        const tensor M(V & D & V.T());
        return symmTensor
        (
            M.xx(), M.xy(), M.xz(),
                    M.yy(), M.yz(),
                            M.zz()
        );
    }

    // Gather a per-rank field into one global field, concatenated in
    // processor order and replicated on every rank.
    template<class T>
    List<T> gatherGlobal(const UList<T>& local)
    {
        List<List<T>> all(Pstream::nProcs());
        all[Pstream::myProcNo()] = local;
        Pstream::gatherList(all);
        Pstream::scatterList(all);
        label n = 0;
        forAll(all, p) n += all[p].size();
        List<T> g(n);
        label idx = 0;
        forAll(all, p) forAll(all[p], i) g[idx++] = all[p][i];
        return g;
    }

    // Bounded, relaxed scale factor (num/den)^expo
    inline scalar boundedScale
    (
        scalar num, scalar den, scalar expo, scalar lo, scalar hi
    )
    {
        const scalar ratio = (mag(den) > SMALL) ? num/den : scalar(1);
        const scalar f = (ratio > 0) ? Foam::pow(ratio, expo) : scalar(1);
        return max(lo, min(hi, f));
    }

    // Area-weighted, component-normalised RMS residual:
    //   sum_c sqrt(sum a*diff_c^2/sum a) / max(max_j|ref_c|, 1).
    // The per-component scale is taken from the reference target.
    inline scalar componentRmsResidual
    (
        const scalarList& area,
        const Field<symmTensor>& diff,
        const Field<symmTensor>& ref
    )
    {
        scalar acc[6] = {0, 0, 0, 0, 0, 0};
        scalar scl[6] = {1, 1, 1, 1, 1, 1};
        scalar A = 0;
        forAll(diff, j)
        {
            A += area[j];
            for (direction d = 0; d < 6; ++d)
            {
                acc[d] += area[j]*sqr(diff[j][d]);
                scl[d] = max(scl[d], mag(ref[j][d]));
            }
        }
        A = max(A, SMALL);
        // Sum normals (xx,yy,zz) then shears (xy,xz,yz).
        const direction ord[6] = {0, 3, 5, 1, 2, 4};
        scalar total = 0;
        for (int k = 0; k < 6; ++k)
        {
            const direction d = ord[k];
            total += Foam::sqrt(acc[d]/A)/scl[d];
        }
        return total;
    }

    // Index of the first element greater than u in a sorted array.
    inline label upperBound(const scalar* l, const label n, const scalar u)
    {
        label lo = 0, hi = n;
        while (hi > lo)
        {
            const label m = (lo + hi)/2;
            if (l[m] <= u) lo = m + 1; else hi = m;
        }
        return lo;
    }

    // Binary search of a cumulative-area array for area-weighted face picking.
    inline label sampleCDF(const scalarList& cumA, scalar u)
    {
        label lo = 0, hi = cumA.size();
        while (hi - lo > 1)
        {
            const label m = (lo + hi)/2;
            if (cumA[m] <= u) lo = m; else hi = m;
        }
        return min(hi, cumA.size() - 1);
    }

}

Foam::tmp<Foam::scalarField>
Foam::turbulentSyntheticEddyInletFvPatchVectorField::rhoOnPatch() const
{
    const surfaceScalarField& phi =
        this->internalField().mesh().lookupObject<surfaceScalarField>("phi");
    if (phi.dimensions() == dimMass/dimTime)
    {
        return tmp<scalarField>
        (
            new scalarField
            (
                patch().lookupPatchField<volScalarField, scalar>("rho")
            )
        );
    }
    return tmp<scalarField>(new scalarField(patch().size(), rhoInf_));
}


void Foam::turbulentSyntheticEddyInletFvPatchVectorField::updateUmean()
{
    const scalarField rho(rhoOnPatch());
    Umean_ = inFlow_->value(db().time().value())
        / (massFlowMode_ ? gSum(rho*patch().magSf()) : 1);
}


Foam::scalar Foam::turbulentSyntheticEddyInletFvPatchVectorField::frictionUtau
(
    scalar Ub,
    scalar nu
) const
{
    // Utau = Ub*sqrt(f/8), f = Darcy friction factor.
    const scalar Re = max(Ub, 1e-4)*delta_/nu;
    scalar sqf;
    if (Re < 3000)
    {
        // Laminar: f = 64/Re, sqf = sqrt(f)
        sqf = sqrt(64.0/Re);
    }
    else
    {
        // Turbulent: Colebrook smooth-pipe, solved by Newton iteration.
        sqf = 0.2;
        for (int i = 0; i < 20; ++i)
        {
            const scalar F  =
                1.93*log10(Re*sqf)*sqf*sqf - sqf - 0.537*sqf*sqf;
            const scalar Fp =
                1.93*sqf/log(10.0)
              + 3.86*sqf*log10(Re*sqf)
              - 1.0 - 1.074*sqf;
            if (mag(Fp) < SMALL) break;
            sqf -= F/Fp;
        }
    }
    return max(Ub, 1e-4)*sqf/sqrt(8.0);
}


// * * * * * * * * * * * * Geometry initialisation * * * * * * * * * * * * * //

void Foam::turbulentSyntheticEddyInletFvPatchVectorField::initializeGeometry()
{
    centerCoord_ = gAverage(patch().Cf());
    const vectorField nf(patch().nf());
    patchNormal_ = -gSum(nf);
    if (mag(patchNormal_) > 0)
    {
        patchNormal_ /= mag(patchNormal_);
    }

    // sqrt(planarError) is the worst-case angular deviation (rad) of a face
    // normal from the patch-averaged normal.
    const scalar planarError =
        returnReduce(max(magSqr(patchNormal_ + nf)), maxOp());
    if (planarError > SMALL)
    {
        const scalar angleDeg =
            Foam::sqrt(planarError)*180.0/constant::mathematical::pi;
        InfoInFunction
            << "Patch " << patch().name() << " is not planar"
            << " (worst face deviates " << angleDeg
            << " deg from patch-averaged normal)" << endl;
    }

    if (!swirlAxisIsInDict_)
    {
        swirlAxis_ = patchNormal_;
    }
    else if (mag(swirlAxis_) > 0)
    {
        swirlAxis_ /= mag(swirlAxis_);
    }

    if (!swirlCenterIsInDict_)
    {
        swirlCenter_ = centerCoord_;
    }

    const polyBoundaryMesh& bm =
        patch().poly().boundaryMesh();

    if (!wallPatchesIsInDict_)
    {
        DynamicList<word> autoNames;
        forAll(bm, patchI)
        {
            if (isA<wallPolyPatch>(bm[patchI]))
            {
                autoNames.append(bm[patchI].name());
            }
        }
        wallPatchNames_ = autoNames;
    }

    computeWallDistance();

    const scalar maxWd = gMax(wallDist_);
    if (maxWd < 0.5*GREAT)
    {
        delta_ = 2*maxVertexWallDist_;
    }
    else
    {
        const vectorField rVec(patch().Cf() - centerCoord_);
        delta_ = 2*gMax(mag(rVec));
    }
    if (delta_ < SMALL)
    {
        delta_ = 1;
    }

    const surfaceScalarField& phi =
        this->internalField().mesh().lookupObject<surfaceScalarField>("phi");
    if (phi.dimensions() == dimMass/dimTime)
    {
        const fvPatchField<scalar>& mu =
            patch().lookupPatchField<volScalarField, scalar>("mu");
        const fvPatchField<scalar>& rho =
            patch().lookupPatchField<volScalarField, scalar>("rho");
        avNu_ = gAverage(scalarField(mu/rho));
    }
    else
    {
        const fvPatchField<scalar>& nu =
            patch().lookupPatchField<volScalarField, scalar>("nu");
        avNu_ = gAverage(nu);
    }

    // Unless specified, compute Utau from the mean inflow.
    if (!UtauIsInDict_ || Utau_ <= 0)
    {
        const scalar Utau_floor = (UtauIsInDict_ && Utau_ <= 0)
            ? mag(Utau_)
            : 0;

        const scalarField rho(rhoOnPatch());
        const scalar t = db().time().value();
        Umean_ = inFlow_->integral(t, t + 10.0)/10.0
            / (massFlowMode_ ? gSum(rho*patch().magSf()) : 1);

        Utau_ = max(frictionUtau(Umean_, avNu_), Utau_floor);
    }

    // Pick the closest bundled DNS dataset by Re_tau.
    const scalar h = 0.5*delta_;
    const scalar ReTau = Utau_*h/avNu_;
    selectedDataset_ = &dnsDatasets[0];
    for (int k = 1; k < dnsDatasetsN; ++k)
    {
        if
        (
            mag(scalar(dnsDatasets[k].Re_tau) - ReTau)
          < mag(scalar(selectedDataset_->Re_tau) - ReTau)
        )
        {
            selectedDataset_ = &dnsDatasets[k];
        }
    }

    // yPlusScale = y+ lookup factor. Floored so DNS profile spans wallDist=h
    const scalar yPlusScalePhys = Utau_/avNu_;
    const scalar yPlusScaleFloor = scalar(selectedDataset_->Re_tau)/h;
    yPlusScale_ = max(yPlusScalePhys, yPlusScaleFloor);

    if (debug)
    {
        Info<< "Patch " << patch().name()
            << ": delta=" << delta_
            << " Utau=" << Utau_
            << " Re_tau=" << ReTau
            << " yPlusScale=" << yPlusScale_
            << " (DNS dataset Re_tau=" << selectedDataset_->Re_tau
            << (yPlusScale_ > yPlusScalePhys + SMALL ? ", yPlusScale floored for full-domain coverage" : "")
            << ")" << endl;
    }
}


void Foam::turbulentSyntheticEddyInletFvPatchVectorField::computeWallDistance()
{
    const polyBoundaryMesh& bm =
        patch().poly().boundaryMesh();

    DynamicList<point> wallCfL;
    DynamicList<vector> wallNfL;
    DynamicList<word> missingNames;
    forAll(wallPatchNames_, k)
    {
        const label patchI = bm.findIndex(wallPatchNames_[k]);
        if (patchI < 0)
        {
            missingNames.append(wallPatchNames_[k]);
            continue;
        }
        const polyPatch& pp = bm[patchI];
        const vectorField cf(pp.faceCentres());
        const vectorField sf(pp.faceAreas());
        forAll(cf, fI)
        {
            wallCfL.append(cf[fI]);
            const scalar mn = mag(sf[fI]);
            wallNfL.append(mn > SMALL ? sf[fI]/mn : sf[fI]);
        }
    }
    if (Pstream::master() && !missingNames.empty() && wallPatchesIsInDict_)
    {
        FatalErrorInFunction
            << "Patch '" << patch().name() << "': wallPatches entries "
            << "not found in the mesh: " << missingNames
            << ". Check spelling against constant/polyMesh/boundary." << endl;
    }

    // Gather wall faces from every rank
    const pointField  wallCf(gatherGlobal<point>(wallCfL));
    const vectorField wallNf(gatherGlobal<vector>(wallNfL));

    wallDist_.setSize(patch().size());
    wallNormal_.setSize(patch().size());

    const vectorField& Cf = patch().Cf();

    if (wallCf.empty())
    {
        if (Pstream::master())
        {
            WarningInFunction
                << "Patch '" << patch().name() << "': no wall faces "
                << "found. " << (wallPatchesIsInDict_
                    ? "User-supplied wallPatches matched nothing."
                    : "Auto-detection found no wallPolyPatch. If your "
                      "walls are declared as type 'patch' you must list them explicitly "
                      "via the `wallPatches` dict entry.")
                << " BC will run without wall scaling. y+ lookups "
                << "saturate at the DNS far-wall value and the near-wall "
                << "Reynolds-stress profile might be wrong." << endl;
        }
        vector vDefault = vector(0, 1, 0) ^ patchNormal_;
        if (mag(vDefault) < 1e-5)
        {
            vDefault = vector(0, 0, 1) ^ patchNormal_;
        }
        vDefault /= mag(vDefault);
        forAll(Cf, i)
        {
            wallDist_[i] = GREAT;
            wallNormal_[i] = vDefault;
        }
        return;
    }

    forAll(Cf, i)
    {
        // Nearest wall face by distance to its centre.
        scalar dminSqr = GREAT;
        label nearestW = -1;
        forAll(wallCf, w)
        {
            const scalar dSqr = magSqr(Cf[i] - wallCf[w]);
            if (dSqr < dminSqr)
            {
                dminSqr = dSqr;
                nearestW = w;
            }
        }

        if (nearestW < 0)
        {
            wallDist_[i] = GREAT;
            wallNormal_[i] = patchNormal_;
            continue;
        }

        // Perpendicular distance to that face's plane.
        const scalar dEuclid = Foam::sqrt(dminSqr);
        const scalar dPerp = max
        (
            scalar(0),
            (Cf[i] - wallCf[nearestW]) & (-wallNf[nearestW])
        );

        // Average inward normals of near-equidistant wall faces (corners).
        const scalar tol = max(1.001*dEuclid, 1e-12);
        vector avgN(Zero);
        forAll(wallCf, w)
        {
            if (mag(Cf[i] - wallCf[w]) <= tol)
            {
                avgN += -wallNf[w];
            }
        }
        if (mag(avgN) < 1e-6)
        {
            avgN = -wallNf[nearestW];
        }

        wallDist_[i] = dPerp;
        wallNormal_[i] = avgN/max(mag(avgN), SMALL);
    }

    // Vertex-based delta
    const pointField& verts = patch().poly().localPoints();
    scalar maxVw = 0;
    forAll(verts, v)
    {
        scalar dminSqr = GREAT;
        label nearestW = -1;
        forAll(wallCf, w)
        {
            const scalar dSqr = magSqr(verts[v] - wallCf[w]);
            if (dSqr < dminSqr) { dminSqr = dSqr; nearestW = w; }
        }
        if (nearestW >= 0)
        {
            const scalar dPerp = max
            (
                scalar(0),
                (verts[v] - wallCf[nearestW]) & (-wallNf[nearestW])
            );
            if (dPerp > maxVw) maxVw = dPerp;
        }
    }
    maxVertexWallDist_ = returnReduce(maxVw, maxOp());
}


// * * * * * * * * * * * * * * * Seed resolution * * * * * * * * * * * * * * //

Foam::label
Foam::turbulentSyntheticEddyInletFvPatchVectorField::resolveSeed
(
    label requested
)
{
    if (requested != -1) return requested;   // fixed seed

    label rnd = 0;
    if (Pstream::master())
    {
        rnd = label(clock::getTime() & 0x7fffffff);
    }
    Pstream::scatter(rnd);   // broadcast master's value to all ranks
    Info<< "Detected seed = -1. Drawing a random seed: " << rnd << endl;
    return rnd;
}


// * * * * * * * * * * * * Per-face eddy template * * * * * * * * * * * * * //

Foam::turbulentSyntheticEddyInletFvPatchVectorField::EddyTemplate
Foam::turbulentSyntheticEddyInletFvPatchVectorField::buildEddyTemplate
(
    const symmTensor& Rlocal,
    scalar Lh,
    const vector& patchNormal,
    const vector& wallNormal,
    scalar Lscale,
    scalar halfWidth,
    scalar cellDx,
    label nCellPerEddy,
    bool useUserR,
    bool useUserL
)
{
    using constant::mathematical::pi;

    EddyTemplate out;

    // Wall-aligned local frame: u = streamwise, v = wall-normal, w = u x v.
    const tensor T(localFrameT(patchNormal, wallNormal));

    // Rotate R to the global frame (DNS path only, user R is already global).
    const symmTensor Rglobal
    (
        useUserR ? Rlocal : localToGlobal(Rlocal, T)
    );

    // Eigendecomposition
    const vector lambda(eigenValues(Rglobal));
    out.Rpg = eigenVectors(Rglobal, lambda).T();
    // Restricting gamma to <= sqrt(8)
    if (lambda[2] <= 0 || lambda[2] > 8.0*(lambda[0] + lambda[1]))
    {
        out.active = false;
        return out;
    }
    out.active = true;
    out.gamma =
        max(scalar(1.0), Foam::sqrt(lambda[2]/(lambda[0] + lambda[1])));
    const scalar g2 = out.gamma*out.gamma;
    out.alphaCoeff = vector
    (
        Foam::sqrt(max(scalar(0), lambda[2] - g2*lambda[0] + g2*lambda[1])),
        Foam::sqrt(max(scalar(0), g2*lambda[0] - g2*lambda[1] + lambda[2])),
        Foam::sqrt(max(scalar(0), g2*lambda[0] + g2*lambda[1] - lambda[2]))
    );

    const scalar Lsc =
        useUserL ? (mag(Lh)*Lscale) : (mag(Lh)*halfWidth*Lscale);
    const scalar localFloor = nCellPerEddy*cellDx;
    out.sigmaX = max(Lsc, localFloor);

    out.c1 = out.gamma/Foam::sqrt(pi*4.0/3.0*pow3(out.sigmaX));

    return out;
}


// * * * * * * * * * * Profile initialisation per face * * * * * * * * * * * //

void Foam::turbulentSyntheticEddyInletFvPatchVectorField::initializeProfiles()
{
    const vectorField& Cf = patch().Cf();
    const scalarField& magSf = patch().magSf();
    const fvMesh& mesh = this->internalField().mesh();
    const label pIdx = patch().index();

    // Lazy init: read custom fields from the start-time directory.
    const word startName = mesh.time().timeName(mesh.time().startTime().value());

    const volSymmTensorField* RUser = nullptr;
    if (!customRName_.empty())
    {
        if (!mesh.foundObject<volSymmTensorField>(customRName_))
        {
            // Not registered yet: read and store from the start-time directory.
            typeIOobject<volSymmTensorField> io
            (
                customRName_, startName, mesh,
                IOobject::MUST_READ, IOobject::NO_WRITE
            );
            if (io.headerOk())
            {
                regIOobject::store(new volSymmTensorField(io, mesh));
            }
        }
        if (mesh.foundObject<volSymmTensorField>(customRName_))
        {
            RUser = &mesh.lookupObject<volSymmTensorField>(customRName_);
        }
        else
        {
            FatalErrorInFunction
                << "customRName '" << customRName_
                << "' set but no volSymmTensorField of that name is registered"
                << " or present in the " << startName << " directory."
                << exit(FatalError);
        }
    }
    const volScalarField* LUser = nullptr;
    if (!customLName_.empty())
    {
        if (!mesh.foundObject<volScalarField>(customLName_))
        {
            typeIOobject<volScalarField> io
            (
                customLName_, startName, mesh,
                IOobject::MUST_READ, IOobject::NO_WRITE
            );
            if (io.headerOk())
            {
                regIOobject::store(new volScalarField(io, mesh));
            }
        }
        if (mesh.foundObject<volScalarField>(customLName_))
        {
            LUser = &mesh.lookupObject<volScalarField>(customLName_);
        }
        else
        {
            FatalErrorInFunction
                << "customLName '" << customLName_
                << "' set but no volScalarField of that name is registered"
                << " or present in the " << startName << " directory."
                << exit(FatalError);
        }
    }
    const volVectorField* UProfileUser = nullptr;
    if (!customUmeanName_.empty())
    {
        if (!mesh.foundObject<volVectorField>(customUmeanName_))
        {
            typeIOobject<volVectorField> io
            (
                customUmeanName_, startName, mesh,
                IOobject::MUST_READ, IOobject::NO_WRITE
            );
            if (io.headerOk())
            {
                regIOobject::store(new volVectorField(io, mesh));
            }
        }
        if (mesh.foundObject<volVectorField>(customUmeanName_))
        {
            UProfileUser = &mesh.lookupObject<volVectorField>(customUmeanName_);
        }
        else
        {
            FatalErrorInFunction
                << "customUmeanName '" << customUmeanName_
                << "' set but no volVectorField of that name is registered"
                << " or present in the " << startName << " directory."
                << exit(FatalError);
        }
    }

    useUserR_ = (RUser != nullptr);
    useUserL_ = (LUser != nullptr);

    U_.setSize(patch().size(), Zero);
    faceRpg_.setSize(patch().size(), tensor::I);
    faceGamma_.setSize(patch().size(), 1.0);
    faceAlphaCoeff_.setSize(patch().size(), Zero);
    faceActive_.setSize(patch().size(), false);
    faceSigmaX_.setSize(patch().size(), 0);
    faceRlocal_.setSize(patch().size(), Zero);
    faceLlocal_.setSize(patch().size(), 0);

    const DNSDataset& dset = *selectedDataset_;
    // Longitudinal two-point integral length L_uu^x.
    const scalar* const Larr = dset.L;

    forAll(Cf, i)
    {
        const scalar r = wallDist_[i]*yPlusScale_;

        // DNS lookup (interpolate Uaxial, Llocal, Rlocal at y+ = r)
        scalar Uaxial = 0;
        scalar Llocal = 0;
        symmTensor Rlocal(Zero);

        const label j = upperBound(dset.yplus, dset.N, r);

        if (j == 0)
        {
            Uaxial = 0;
            Llocal = Larr[0];
            Rlocal = symmTensor
            (
                dset.Ruu[0], dset.Ruv[0], dset.Ruw[0],
                             dset.Rvv[0], dset.Rvw[0],
                                          dset.Rww[0]
            );
        }
        else if (j >= dset.N)
        {
            Uaxial = dset.U[dset.N - 1];
            Llocal = Larr[dset.N - 1];
            Rlocal = symmTensor
            (
                dset.Ruu[dset.N - 1], dset.Ruv[dset.N - 1], dset.Ruw[dset.N - 1],
                                      dset.Rvv[dset.N - 1], dset.Rvw[dset.N - 1],
                                                            dset.Rww[dset.N - 1]
            );
        }
        else
        {
            const scalar dr = dset.yplus[j] - dset.yplus[j - 1];
            const scalar w0 = (dset.yplus[j] - r)/dr;
            const scalar w1 = (r - dset.yplus[j - 1])/dr;
            Uaxial = dset.U[j - 1]*w0 + dset.U[j]*w1;
            Llocal = Larr[j - 1]*w0 + Larr[j]*w1;
            Rlocal =
                symmTensor
                (
                    dset.Ruu[j - 1], dset.Ruv[j - 1], dset.Ruw[j - 1],
                                     dset.Rvv[j - 1], dset.Rvw[j - 1],
                                                      dset.Rww[j - 1]
                )*w0
              + symmTensor
                (
                    dset.Ruu[j], dset.Ruv[j], dset.Ruw[j],
                                 dset.Rvv[j], dset.Rvw[j],
                                              dset.Rww[j]
                )*w1;
        }

        // U setup
        if (UProfileUser)
        {
            U_[i] = UProfileUser->boundaryField()[pIdx][i];
        }
        else
        {
            U_[i] = Uaxial*patchNormal_;
        }

        // DNS R/L: wall-aligned local frame, dimensionless. User R/L: taken
        // as-is (global frame, physical units).
        faceRlocal_[i] = RUser ? RUser->boundaryField()[pIdx][i] : Rlocal;
        faceLlocal_[i] = LUser ? LUser->boundaryField()[pIdx][i] : Llocal;
    }

    rebuildTemplates();
    v0_ = 2*gSum(magSf)*maxSigmaX_;

    // Snapshot the reference target before the MC overwrites faceRlocal_.
    // The downstream loop drives the realized R toward this.
    if (ds_.distance > 0) faceRlocalRef_ = faceRlocal_;

    applyInletRcorrectionMC();  // MC iteration on faceRlocal_ + re-runs

    scalar Uavg = gSum(U_ & -patch().Sf())/gSum(magSf);
    if (mag(Uavg) < SMALL)
    {
        Uavg = 1;
    }
    U_ *= 1.0/Uavg;
    updateUmean();

    if (mag(swirlNumber_) > SMALL)
    {
        const scalarField rho(rhoOnPatch());
        scalarField rPerp(patch().size(), 0);
        vectorField rPerpVec(patch().size(), Zero);
        scalarField UaxialF(patch().size(), 0);

        forAll(Cf, i)
        {
            const vector d = Cf[i] - swirlCenter_;
            rPerpVec[i] = d - (d & swirlAxis_)*swirlAxis_;
            rPerp[i] = mag(rPerpVec[i]);
            UaxialF[i] = U_[i] & swirlAxis_;
        }

        const scalar Rm = gSum(rPerp*magSf)/gSum(magSf);
        const scalar G1 =
            gSum(rho*sqr(UaxialF)*magSf*sqr(rPerp));
        const scalar G2 = gSum(rho*sqr(UaxialF)*magSf);

        if (G1 > SMALL && G2 > SMALL)
        {
            const scalar swirlFrac = swirlNumber_*Rm*G2/G1;
            forAll(Cf, i)
            {
                vector tang(Zero);
                if (rPerp[i] > SMALL)
                {
                    tang = (swirlAxis_ ^ rPerpVec[i])/rPerp[i];
                }
                U_[i] += UaxialF[i]*swirlFrac*rPerp[i]*tang;
            }
        }
    }
}


// * * * * * * * * * Eddy-template construction  * * * * * * * * * * //

void Foam::turbulentSyntheticEddyInletFvPatchVectorField::rebuildTemplates()
{
    // Per-face buildEddyTemplate from the current faceRlocal_ / faceLlocal_.
    // Sets faceRpg_, faceActive_, faceGamma_, faceAlphaCoeff_, faceSigmaX_,
    // plus the global maxSigmaX_ summary.
    const vectorField& Cf = patch().Cf();
    const scalarField& magSf = patch().magSf();
    // Cap cellDx at half the inlet size (guards degenerate faces).
    const scalar cellDxCap = 0.5*Foam::sqrt(max(gSum(magSf), SMALL));
    const scalarField cellDx
    (
        min
        (
            max(Foam::sqrt(magSf), 2/patch().deltaCoeffs()),
            cellDxCap
        )
    );
    const scalar h = 0.5*delta_;

    forAll(Cf, i)
    {
        const EddyTemplate et = buildEddyTemplate
        (
            faceRlocal_[i], faceLlocal_[i],
            patchNormal_, wallNormal_[i],
            Lscale_, h, cellDx[i],
            nCellPerEddy_,
            useUserR_, useUserL_
        );
        faceRpg_[i] = et.Rpg;
        faceActive_[i] = et.active;
        faceGamma_[i] = et.gamma;
        faceAlphaCoeff_[i] = et.alphaCoeff;
        faceSigmaX_[i] = et.sigmaX;
    }
    maxSigmaX_ = gMax(faceSigmaX_);
}


// * * * * * * * * * Shared R-correction sample set * * * * * * * * * * * * * //

void Foam::turbulentSyntheticEddyInletFvPatchVectorField::buildSampleSet
(
    RcSampleSet& ss
) const
{
    ss.valid = false;
    const scalar sampleAlpha = 0.2;   // Poisson-disk r_excl = alpha*sigma_x(face)
    // Near-wall exclusion: drop faces within wnExcludeCoeff*(sigma_x/gamma) of
    // the wall, where the eddy is larger than its wall distance and the
    // realized R cannot follow the target's decline to zero.
    const scalar wnExcludeCoeff = ss.wnExcludeCoeff;

    // ----- 1. Gather active face state onto every rank -----
    const vectorField& patchCf = patch().Cf();
    const scalarField& patchMagSf = patch().magSf();
    DynamicList<vector> CfLocal, wnLocal;
    DynamicList<scalar> wdLocal, aLocal, LLocal;
    DynamicList<symmTensor> RLocal;
    DynamicList<label> origIdxLocal;
    forAll(wallDist_, i)
    {
        if (!faceActive_[i]) continue;
        CfLocal.append(patchCf[i]);
        wdLocal.append(wallDist_[i]);
        aLocal.append(patchMagSf[i]);
        wnLocal.append(wallNormal_[i]);
        RLocal.append(faceRlocal_[i]);
        LLocal.append(faceLlocal_[i]);
        origIdxLocal.append(i);
    }
    const List<vector> CfF(gatherGlobal<vector>(CfLocal));
    const List<vector> wnF(gatherGlobal<vector>(wnLocal));
    const scalarList wdF(gatherGlobal<scalar>(wdLocal));
    const scalarList aF(gatherGlobal<scalar>(aLocal));
    const scalarList LF(gatherGlobal<scalar>(LLocal));
    const List<symmTensor> RF(gatherGlobal<symmTensor>(RLocal));
    const label nF = wdF.size();
    if (nF < 4) return;
    const scalar A_total = sum(scalarField(aF));

    // This rank's offset into the gathered (proc-order) arrays, for write-back.
    labelList procCounts(Pstream::nProcs(), 0);
    procCounts[Pstream::myProcNo()] = origIdxLocal.size();
    Pstream::gatherList(procCounts);
    Pstream::scatterList(procCounts);
    label off = 0;
    for (label p = 0; p < Pstream::myProcNo(); ++p) off += procCounts[p];

    // ----- 2. Per-face geometry, sigmaX_face, R_ref_face, Tframe_face -----
    const scalar h = 0.5*delta_;
    const scalar cellDxCapMC =
        0.5*Foam::sqrt(max(gSum(patch().magSf()), SMALL));
    const scalarField cellDxLocal
    (
        min
        (
            max(Foam::sqrt(patch().magSf()), 2/patch().deltaCoeffs()),
            cellDxCapMC
        )
    );
    const scalar cellDxAvg = gAverage(cellDxLocal);
    Field<symmTensor> R_ref_face(nF);
    Field<tensor> Tframe_face(nF);
    scalarField sigmaX_face(nF, 0);
    scalarField sigmaWN_face(nF, 0);
    scalar maxSigma = 0;
    for (label i = 0; i < nF; ++i)
    {
        R_ref_face[i] = RF[i];
        Tframe_face[i] = localFrameT(patchNormal_, wnF[i]);
        const EddyTemplate et = buildEddyTemplate
        (
            RF[i], LF[i], patchNormal_, wnF[i],
            Lscale_, h, cellDxAvg,
            nCellPerEddy_,
            useUserR_, useUserL_
        );
        if (et.active)
        {
            sigmaX_face[i] = et.sigmaX;
            sigmaWN_face[i] = et.sigmaX/max(et.gamma, SMALL);
            maxSigma = max(maxSigma, et.sigmaX);
        }
    }
    if (maxSigma <= 0) return;
    if (debug)
    {
        scalar swMin = GREAT, swMax = 0;
        for (label i = 0; i < nF; ++i)
            if (sigmaWN_face[i] > 0)
            { swMin = min(swMin, sigmaWN_face[i]); swMax = max(swMax, sigmaWN_face[i]); }
        Info<< "  near-wall exclude band = " << wnExcludeCoeff << "*(sigmaX/gamma);"
            << " sigmaX/gamma in [" << swMin << ',' << swMax << "]"
            << "  -> excludes wallDist up to ~" << wnExcludeCoeff*swMin
            << " (min) .. " << wnExcludeCoeff*swMax << " (max)" << endl;
    }

    // ----- 3. Poisson-disk sample selection (density proportional to 1/sigma_x^2) -----
    labelList sampleIdx;
    {
        labelList shuffled(nF);
        forAll(shuffled, i) shuffled[i] = i;
        randomGenerator rngPick(randomGenerator::seed(seed_ + 7919));
        for (label i = nF - 1; i > 0; --i)
        {
            Swap(shuffled[i], shuffled[rngPick.sampleAB<label>(0, i + 1)]);
        }
        DynamicList<label> picks;
        forAll(shuffled, oi)
        {
            const label k = shuffled[oi];
            if (sigmaX_face[k] <= 0) continue;
            if (wdF[k] <= wnExcludeCoeff*sigmaWN_face[k]) continue;
            const scalar rHere = sampleAlpha*sigmaX_face[k];
            bool ok = true;
            forAll(picks, p)
            {
                const label q = picks[p];
                const scalar rT = max(rHere, sampleAlpha*sigmaX_face[q]);
                if (magSqr(CfF[k] - CfF[q]) < sqr(rT)) { ok = false; break; }
            }
            if (ok) picks.append(k);
        }
        sampleIdx.transfer(picks);
    }
    const label nS = sampleIdx.size();
    if (nS < 4) return;

    // Sample-level slice arrays.
    List<vector> Cf_s(nS);
    scalarList a_s(nS);
    Field<symmTensor> R_ref_s(nS);
    Field<tensor> Tframe_s(nS);
    forAll(sampleIdx, j)
    {
        const label k = sampleIdx[j];
        Cf_s[j] = CfF[k];
        a_s[j] = aF[k];
        R_ref_s[j] = RF[k];
        Tframe_s[j] = Tframe_face[k];
    }

    // ----- 4. Sparse 2D kernel matrix W (samples x samples) + face<-sample -----
    const scalar mean_face_spacing = Foam::sqrt(A_total/scalar(nF));
    const scalar sigmaKernel = max(maxSigma/1.5, 2.0*mean_face_spacing);

    List<labelList> nbrs_ss(nS);
    List<scalarField> W_ss(nS);
    scalarField wnorm_ss(nS, 0);
    for (label j = 0; j < nS; ++j)
    {
        DynamicList<label> ids;
        DynamicList<scalar> ws;
        for (label k = 0; k < nS; ++k)
        {
            const scalar d = mag(Cf_s[k] - Cf_s[j]);
            if (d >= sigmaKernel) continue;
            const scalar a = d/sigmaKernel;
            const scalar w = pow3(1.0 - a*a);
            if (w <= 0) continue;
            ids.append(k);
            ws.append(w);
        }
        nbrs_ss[j].transfer(ids);
        W_ss[j].transfer(ws);
        scalar wn = 0;
        forAll(nbrs_ss[j], kk) wn += W_ss[j][kk]*a_s[nbrs_ss[j][kk]];
        wnorm_ss[j] = max(wn, SMALL);
    }

    List<labelList> nbrs_fs(nF);
    List<scalarField> W_fs(nF);
    scalarField wnorm_fs(nF, 0);
    for (label i = 0; i < nF; ++i)
    {
        DynamicList<label> ids;
        DynamicList<scalar> ws;
        for (label j = 0; j < nS; ++j)
        {
            const scalar d = mag(Cf_s[j] - CfF[i]);
            if (d >= sigmaKernel) continue;
            const scalar a = d/sigmaKernel;
            const scalar w = pow3(1.0 - a*a);
            if (w <= 0) continue;
            ids.append(j);
            ws.append(w);
        }
        if (ids.empty())   // fallback: nearest sample
        {
            scalar dmin = GREAT; label kmin = 0;
            for (label j = 0; j < nS; ++j)
            {
                const scalar d = mag(Cf_s[j] - CfF[i]);
                if (d < dmin) { dmin = d; kmin = j; }
            }
            ids.append(kmin); ws.append(1.0);
        }
        nbrs_fs[i].transfer(ids);
        W_fs[i].transfer(ws);
        scalar wn = 0;
        forAll(W_fs[i], kk) wn += W_fs[i][kk];
        wnorm_fs[i] = max(wn, SMALL);
    }

    // ----- store into the sample set -----
    ss.nF = nF; ss.A_total = A_total;
    ss.CfF = CfF; ss.wnF = wnF; ss.wdF = wdF; ss.aF = aF; ss.LF = LF;
    ss.h = h; ss.cellDxAvg = cellDxAvg; ss.maxSigma = maxSigma;
    ss.R_ref_face = R_ref_face; ss.Tframe_face = Tframe_face;
    ss.sigmaX_face = sigmaX_face;
    ss.nS = nS;
    ss.Cf_s = Cf_s;
    ss.a_s = a_s;
    ss.R_ref_s = R_ref_s; ss.Tframe_s = Tframe_s;
    ss.sigmaKernel = sigmaKernel;
    ss.nbrs_ss = nbrs_ss; ss.W_ss = W_ss; ss.wnorm_ss = wnorm_ss;
    ss.nbrs_fs = nbrs_fs; ss.W_fs = W_fs; ss.wnorm_fs = wnorm_fs;
    ss.origIdxLocal = origIdxLocal; ss.off = off;
    ss.valid = true;
}


void Foam::turbulentSyntheticEddyInletFvPatchVectorField::scatterSampleToFace
(
    const RcSampleSet& ss,
    const Field<symmTensor>& Rs,
    Field<symmTensor>& Rface
) const
{
    for (label i = 0; i < ss.nF; ++i)
    {
        symmTensor acc = symmTensor::zero;
        forAll(ss.nbrs_fs[i], kk)
            acc += ss.W_fs[i][kk]*Rs[ss.nbrs_fs[i][kk]];
        Rface[i] = acc/ss.wnorm_fs[i];
    }
}


Foam::scalar Foam::turbulentSyntheticEddyInletFvPatchVectorField::rmsResidual
(
    const RcSampleSet& ss,
    const Field<symmTensor>& diff
) const
{
    return componentRmsResidual(ss.a_s, diff, ss.R_ref_s);
}


// * * * * * * * inlet MC R-correction helpers * * * * * * * * * * //

Foam::vector
Foam::turbulentSyntheticEddyInletFvPatchVectorField::eddyImageContribution
(
    const syntheticEddy& ed, const point& pos, scalar reach, const point& xp
) const
{
    vector up(Zero);
    if (magSqr(xp - pos) <= sqr(reach)) up += ed.uPrime(xp, patchNormal_);
    forAll(cyclicTranslations_, c)
    {
        const vector& t = cyclicTranslations_[c];
        vector sh = xp + t;
        if (magSqr(sh - pos) <= sqr(reach)) up += ed.uPrime(sh, patchNormal_);
        sh = xp - t;
        if (magSqr(sh - pos) <= sqr(reach)) up += ed.uPrime(sh, patchNormal_);
    }
    return up;
}


void
Foam::turbulentSyntheticEddyInletFvPatchVectorField::spawnEddyBatch
(
    const RcSampleSet& ss,
    const scalarList& cumA,
    randomGenerator& rng,
    label nEddies,
    const Field<symmTensor>& R_target_face,
    symmTensorField& Rsum
) const
{
    for (label n = 0; n < nEddies; ++n)
    {
        const label fi = sampleCDF(cumA, rng.scalarAB(0, cumA.last()));
        const symmTensor RlocalAtEddy = R_target_face[fi];
        const EddyTemplate et = buildEddyTemplate
        (
            RlocalAtEddy, ss.LF[fi], patchNormal_, ss.wnF[fi],
            Lscale_, ss.h, ss.cellDxAvg,
            nCellPerEddy_,
            useUserR_, useUserL_
        );
        if (!et.active) continue;
        const vector sigmaE
        (
            et.sigmaX/et.gamma, et.sigmaX/et.gamma, et.sigmaX
        );
        const vector alpha
        (
            (rng.scalar01() < 0.5 ? 1.0 : -1.0)*et.alphaCoeff[0]/et.sigmaX,
            (rng.scalar01() < 0.5 ? 1.0 : -1.0)*et.alphaCoeff[1]/et.sigmaX,
            (rng.scalar01() < 0.5 ? 1.0 : -1.0)*et.alphaCoeff[2]/et.sigmaX
        );
        const point pos
        (
            rng.scalarAB(-ss.maxSigma, ss.maxSigma),
            ss.CfF[fi].y(),
            ss.CfF[fi].z()
        );
        syntheticEddy ed(pos, 0, sigmaE, alpha, et.Rpg, et.c1, scalar(0));
        const scalar reach = et.sigmaX;

        // Remove each eddy's area-weighted patch mean before accumulating
        vector me(Zero);
        for (label fm = 0; fm < ss.nF; ++fm)
            me += ss.aF[fm]*eddyImageContribution(ed, pos, reach, ss.CfF[fm]);
        me /= max(ss.A_total, SMALL);

        for (label jj = 0; jj < ss.nS; ++jj)
            Rsum[jj] += sqr(eddyImageContribution(ed, pos, reach, ss.Cf_s[jj]) - me);
    }
}


void
Foam::turbulentSyntheticEddyInletFvPatchVectorField::writeRcDebugField
(
    const word& name,
    const Field<symmTensor>& faceLocalDict,
    const RcSampleSet& ss,
    scalar physScale
) const
{
    const fvMesh& mesh = this->internalField().mesh();
    const label pIdx = patch().index();
    const word t0 = mesh.time().timeName(mesh.time().startTime().value());

    volSymmTensorField fld
    (
        IOobject
        (
            name, t0, mesh,
            IOobject::NO_READ, IOobject::NO_WRITE, false
        ),
        mesh,
        dimensionedSymmTensor(sqr(dimVelocity), Zero)
    );
    forAll(ss.origIdxLocal, k)
    {
        const label gf = ss.off + k;          // global gathered index
        const label pf = ss.origIdxLocal[k];   // patch-local face index
        fld.boundaryFieldRef()[pIdx][pf] =
            localToGlobal(faceLocalDict[gf], ss.Tframe_face[gf])*physScale;
    }
    // Force the startTime directory (bypass updateInstance).
    fld.instance() = t0;
    fileHandler().writeObject
    (
        fld, mesh.time().writeFormat(), IOstream::currentVersion,
        mesh.time().writeCompression()
    );
}


// * * * * * * * inlet MC R-correction * * * * * * * * * * * * * * //

void Foam::turbulentSyntheticEddyInletFvPatchVectorField::applyInletRcorrectionMC()
{
    if (!inletRcorrection_) return;

    const clockTime timer;

    // Internal tunables
    const scalar kAdd        = 0.3;   // Landweber step size
    const label  N_max       = 5000000;
    const scalar gammaCap    = 3.0;
    const label  nIter_max   = 30;
    const label  rollK       = 5;
    const label  patience    = 3;

    // ----- 1. Gather, Poisson-disk sampling and kernel maps -----
    RcSampleSet ss;
    buildSampleSet(ss);
    if (!ss.valid) return;

    // Read-only aliases into the sample set. Fields consumed only inside
    // spawnEddyBatch/writeRcDebugField are reached via ss there, not aliased here.
    const label nF = ss.nF;
    const scalar A_total = ss.A_total;
    const scalarList& aF = ss.aF;
    const scalar maxSigma = ss.maxSigma;
    const Field<symmTensor>& R_ref_face = ss.R_ref_face;
    const scalarField& sigmaX_face = ss.sigmaX_face;
    const label nS = ss.nS;
    const scalarList& a_s = ss.a_s;
    const Field<symmTensor>& R_ref_s = ss.R_ref_s;
    const Field<tensor>& Tframe_s = ss.Tframe_s;
    const scalar sigmaKernel = ss.sigmaKernel;
    const List<labelList>& nbrs_ss = ss.nbrs_ss;
    const List<scalarField>& W_ss = ss.W_ss;
    const scalarField& wnorm_ss = ss.wnorm_ss;
    const labelList& origIdxLocal = ss.origIdxLocal;
    const label off = ss.off;

    // ----- 2. Auto-N + per-rank eddy counts -----
    // N from the area-weighted mean eddy size, N = avgOverlaps*V_box/V_eddy_avg.
    const scalar V_box = 2.0*maxSigma*A_total;
    scalar S3_A = 0, A_active = 0;
    scalar sigmaMin = GREAT;
    label nActive = 0;
    forAll(sigmaX_face, i)
    {
        if (sigmaX_face[i] > 0)
        {
            S3_A += pow3(sigmaX_face[i])*aF[i];
            A_active += aF[i];
            sigmaMin = min(sigmaMin, sigmaX_face[i]);
            ++nActive;
        }
    }
    if (nActive == 0 || A_active <= 0) return;
    const scalar sigmaCubeAvg = S3_A/A_active;
    const scalar Veddy_avg    = (4.0/3.0)*pi*sigmaCubeAvg/sqr(gammaCap);
    const scalar avgOverlaps  = mcOverlaps_;  // eddy-overlap density (dict knob)
    label N = label(min(scalar(N_max),
                        max(scalar(1000),
                            avgOverlaps*V_box/max(Veddy_avg, SMALL))));

    scalarList cumA(nF, scalar(0));
    { scalar c = 0; forAll(cumA, k) { c += aF[k]; cumA[k] = c; } }

    // Rank-specific seed: independent eddy streams, reduced each iter.
    randomGenerator rng
    (
        randomGenerator::seed(seed_ + Pstream::myProcNo()*7919)
    );
    // Per-rank eddy count: rank 0 takes the remainder so sum N_local = N.
    const label nProcs_ = Pstream::nProcs();
    const label N_local =
        N/nProcs_ + (Pstream::myProcNo() == 0 ? N % nProcs_ : 0);

    // The spawn already produces dimensionless R. physScale only renders the
    // optional debug fields in physical units.
    const scalar physScale = useUserR_ ? Rscale_ : sqr(Utau_)*Rscale_;

    if (debug)
    {
        Info<< "inletRcorrection (MC) on patch '" << patch().name()
            << "': N=" << N;
        if (Pstream::parRun())
            Info<< " (N_local~" << N/nProcs_ << "x" << nProcs_ << " ranks)";
        Info<< "  nF=" << nF << "  nS=" << nS
            << "  sigmaKernel=" << sigmaKernel
            << "  sigmaX(min,<>_A^1/3,max)=("
            << sigmaMin << "," << Foam::cbrt(sigmaCubeAvg) << ","
            << maxSigma << ")" << endl;
    }

    // c2 normalises sum N spawned eddies. N is fixed, so c2 is constant.
    const scalar c2 = Foam::sqrt((315.0/16.0)*V_box/scalar(N));

    // ----- 3. MC iteration loop -----
    Field<symmTensor> R_target_s = R_ref_s;
    Field<symmTensor> R_target_face = R_ref_face;
    Field<symmTensor> R_target_s_best = R_target_s;
    // Realized R from the initial reference target, for the R_init_realized debug field.
    Field<symmTensor> Rrl_face_init(nF, Zero);
    scalarField resRms_buf(rollK, GREAT);
    label rb_idx = 0;
    scalar bestRms = GREAT;
    label bestIter = 0;
    label noImprove = 0;

    for (label iter = 0; iter < nIter_max; ++iter)
    {
        symmTensorField Rsum(nS, Zero);
        spawnEddyBatch(ss, cumA, rng, N_local, R_target_face, Rsum);
        // Allreduce across ranks (sum N_local = N, so c2 normalises correctly).
        if (Pstream::parRun())
        {
            forAll(Rsum, j) reduce(Rsum[j], sumOp());
        }

        Field<symmTensor> Rrl_s(nS), diff_s(nS);
        for (label j = 0; j < nS; ++j)
        {
            Rrl_s[j] = globalToLocal(c2*c2*Rsum[j], Tframe_s[j]);
            diff_s[j] = Rrl_s[j] - R_ref_s[j];
        }
        const scalar resRms = rmsResidual(ss, diff_s);

        // Landweber step on the kernel-smoothed residual, then re-impose
        // realizability.
        {
            Field<symmTensor> resid_s(nS);
            for (label j = 0; j < nS; ++j)
                resid_s[j] = R_ref_s[j] - Rrl_s[j];
            for (label j = 0; j < nS; ++j)
            {
                symmTensor acc = symmTensor::zero;
                forAll(nbrs_ss[j], kk)
                {
                    const label k = nbrs_ss[j][kk];
                    acc += W_ss[j][kk]*a_s[k]*resid_s[k];
                }
                R_target_s[j] = R_target_s[j] + kAdd*(acc/wnorm_ss[j]);
            }
            for (label j = 0; j < nS; ++j)
                R_target_s[j] = projectPDjac(R_target_s[j]);
        }

        // Scatter sample -> face for next iter's spawn.
        scatterSampleToFace(ss, R_target_s, R_target_face);
        if (debug && iter == 0)
        {
            scatterSampleToFace(ss, Rrl_s, Rrl_face_init);
        }

        if (debug) Info<< "  iter " << iter << "  RMS=" << resRms << endl;

        if (resRms < bestRms)
        {
            bestRms = resRms;
            bestIter = iter;
            R_target_s_best = R_target_s;
        }
        resRms_buf[rb_idx] = resRms;
        rb_idx = (rb_idx + 1) % rollK;
        if (iter >= rollK && iter - bestIter >= patience)
        {
            noImprove++;
            if (noImprove >= patience)
            {
                if (debug)
                    Info<< "  early stop at iter " << iter
                        << " (best RMS=" << bestRms << " at iter "
                        << bestIter << ")" << endl;
                break;
            }
        }
        else { noImprove = 0; }
    }

    // ----- 4. Scatter the best-fit target back to every local face -----
    Field<symmTensor> R_target_face_best(nF);
    scatterSampleToFace(ss, R_target_s_best, R_target_face_best);
    forAll(origIdxLocal, k)
    {
        faceRlocal_[origIdxLocal[k]] = R_target_face_best[off + k];
    }

    // ----- 5. Debug: write initial/final target and realized R as
    // volSymmTensorFields (global frame, physical units), zero off the patch.
    if (debug)
    {
        symmTensorField Rsum(nS, Zero);
        spawnEddyBatch(ss, cumA, rng, N_local, R_target_face_best, Rsum);
        if (Pstream::parRun())
        {
            forAll(Rsum, j) reduce(Rsum[j], sumOp());
        }
        Field<symmTensor> Rrl_s_scaled(nS), diff_scaled(nS);
        for (label j = 0; j < nS; ++j)
        {
            Rrl_s_scaled[j] = globalToLocal(c2*c2*Rsum[j], Tframe_s[j]);
            diff_scaled[j] = Rrl_s_scaled[j] - R_ref_s[j];
        }
        Info<< "  final realized: RMS=" << rmsResidual(ss, diff_scaled) << endl;

        Field<symmTensor> Rrl_face_scaled(nF);
        scatterSampleToFace(ss, Rrl_s_scaled, Rrl_face_scaled);

        // Nonzero only on this patch: rotate the dict-frame tensor to global
        // and scale to physical.
        writeRcDebugField("R_init_target",    R_ref_face,         ss, physScale);
        writeRcDebugField("R_init_realized",  Rrl_face_init,      ss, physScale);
        writeRcDebugField("R_final_target",   R_target_face_best, ss, physScale);
        writeRcDebugField("R_final_realized", Rrl_face_scaled,    ss, physScale);
    }

    rebuildTemplates();
    if (debug)
        Info<< "inletRcorrection (MC) done.  Wall time: "
            << timer.elapsedTime()
            << " s  (N=" << N << " nS=" << nS << " bestIter=" << bestIter
            << " bestRms=" << bestRms << ")" << endl;
}


// * * * * * * Downstream-measurement R-correction (runtime loop) * * * * * * //

void
Foam::turbulentSyntheticEddyInletFvPatchVectorField::setupDownstreamProbes()
{
    dsState_.setupDone = true;
    if (ds_.distance <= 0) return;

    // ----- 1. Gather every active face onto every rank -----
    const vectorField& patchCf = patch().Cf();
    const scalarField& patchMagSf = patch().magSf();
    const bool haveRef = (faceRlocalRef_.size() == faceRlocal_.size());
    DynamicList<vector> CfL, wnL;
    DynamicList<scalar> wdL, aL, sxL;
    DynamicList<symmTensor> RL, RrefL;
    forAll(wallDist_, i)
    {
        if (!faceActive_[i]) continue;
        CfL.append(patchCf[i]);
        wdL.append(wallDist_[i]);
        aL.append(patchMagSf[i]);
        wnL.append(wallNormal_[i]);
        sxL.append(faceSigmaX_[i]);
        RL.append(faceRlocal_[i]);
        RrefL.append(haveRef ? faceRlocalRef_[i] : faceRlocal_[i]);
    }
    const List<vector> CfF(gatherGlobal<vector>(CfL));
    const List<vector> wnF(gatherGlobal<vector>(wnL));
    const scalarList wdF(gatherGlobal<scalar>(wdL));
    const scalarList aF(gatherGlobal<scalar>(aL));
    const scalarList sxF(gatherGlobal<scalar>(sxL));
    const List<symmTensor> RF(gatherGlobal<symmTensor>(RL));
    const List<symmTensor> RrefF(gatherGlobal<symmTensor>(RrefL));
    const label nP = wdF.size();
    if (nP < 4)
    {
        WarningInFunction
            << "downstream R-correction: too few active faces on patch '"
            << patch().name() << "'; disabling." << endl;
        ds_.distance = 0;
        return;
    }

    // ----- 2. Probes + parallel ownership -----
    dsState_.probeX.setSize(nP);
    dsState_.probeWd.setSize(nP);
    dsState_.probeArea.setSize(nP);
    dsState_.probeWn.setSize(nP);
    dsState_.maxSigma = 0;
    for (label p = 0; p < nP; ++p)
    {
        dsState_.probeX[p] = CfF[p] + ds_.distance*patchNormal_;
        dsState_.probeWd[p] = wdF[p];
        dsState_.probeArea[p] = aF[p];
        dsState_.probeWn[p] = wnF[p];
        dsState_.maxSigma = max(dsState_.maxSigma, sxF[p]);
    }

    // Lowest rank containing a probe owns it.
    // Probes nobody contains are outside of the domain.
    const fvMesh& mesh = this->internalField().mesh();
    const meshSearch& ms = meshSearch::New(mesh);
    dsState_.localCell = labelList(nP, -1);
    dsState_.ownerProc = labelList(nP, labelMax);
    for (label p = 0; p < nP; ++p)
    {
        const label celli = ms.findCell(dsState_.probeX[p]);
        if (celli >= 0)
        {
            dsState_.localCell[p] = celli;
            dsState_.ownerProc[p] = Pstream::myProcNo();
        }
    }
    if (Pstream::parRun())
        forAll(dsState_.ownerProc, p) reduce(dsState_.ownerProc[p], minOp());

    DynamicList<label> mine;
    label nDead = 0;
    for (label p = 0; p < nP; ++p)
    {
        if (dsState_.ownerProc[p] == labelMax)             // nobody found it
        {
            dsState_.localCell[p] = -1;
            ++nDead;
        }
        else if (dsState_.ownerProc[p] != Pstream::myProcNo())
        {
            dsState_.localCell[p] = -1;                    // demote halo duplicate
        }
        else
        {
            mine.append(p);
        }
    }
    dsState_.myProbes.transfer(mine);
    if (nDead > nP/2)
    {
        FatalErrorInFunction
            << "downstream R-correction on patch '" << patch().name()
            << "': " << nDead << '/' << nP << " probe points at Delta="
            << ds_.distance << " fall outside the domain. Reduce "
            << "correctionDistance." << exit(FatalError);
    }

    // ----- 3. Reference profile and its area-weighted mean (local frame) -----
    dsState_.probeRef = Field<symmTensor>(RrefF);
    dsState_.probeTarget = Field<symmTensor>(RF);   // applied target (diagnostics)
    symmTensor refSum(symmTensor::zero), tgtSum(symmTensor::zero);
    scalar aSum = 0;
    forAll(RrefF, p)
    {
        const scalar a = aF[p];
        refSum += a*RrefF[p];
        tgtSum += a*RF[p];
        aSum += a;
    }
    aSum = max(aSum, SMALL);
    dsState_.refMean = refSum/aSum;
    dsState_.appliedMean = tgtSum/aSum;

    // ----- 4. Accumulators + controller state -----
    dsState_.sumU = vectorField(nP, Zero);
    dsState_.sumUU = symmTensorField(nP, symmTensor::zero);
    dsState_.sampleCount = 0;
    dsState_.windowEddies = ds_.windowEddies;   // grows x1.5 on a stall
    dsState_.pass = 0;
    dsState_.ph = downstreamState::phase::warmup;
    dsState_.phaseStartTime = db().time().value();
    dsState_.bestTkeErr = GREAT;
    dsState_.noImprove = 0;
    dsState_.belowTol = 0;
    dsState_.bestTarget = faceRlocal_;

    if (debug)
    {
        label nMine = dsState_.myProbes.size();
        reduce(nMine, sumOp());
        Info<< "downstream R-correction (uniform): patch '" << patch().name()
            << "' nP=" << nP << " owned=" << nMine << " dead=" << nDead
            << "  correctionDistance=" << ds_.distance << " (="
            << ds_.distance/max(0.5*delta_, SMALL) << " half-lengths)" << endl;
    }
}


void
Foam::turbulentSyntheticEddyInletFvPatchVectorField::stepDownstreamController()
{
    if (ds_.distance <= 0) return;
    if (!dsState_.setupDone) setupDownstreamProbes();
    if (ds_.distance <= 0 || dsState_.ph == downstreamState::phase::done) return;

    const scalar t = db().time().value();
    if (t < ds_.startTime)
    {
        dsState_.phaseStartTime = t;   // measure the warmup from the gate opening
        return;
    }

    const scalar Ub = max(mag(Umean_), SMALL);

    if (dsState_.ph == downstreamState::phase::warmup)
    {
        const scalar warm = 1.5*ds_.distance/Ub;
        if (t - dsState_.phaseStartTime >= warm)
        {
            dsState_.sumU = vector::zero;
            dsState_.sumUU = symmTensor::zero;
            dsState_.sampleCount = 0;
            dsState_.ph = downstreamState::phase::accumulate;
            dsState_.phaseStartTime = t;
        }
        return;
    }

    // accumulate: read U at each owned probe cell.
    const Field<vector>& Uint = this->internalField();
    forAll(dsState_.myProbes, pp)
    {
        const label p = dsState_.myProbes[pp];
        const vector& u = Uint[dsState_.localCell[p]];
        dsState_.sumU[p] += u;
        dsState_.sumUU[p] += sqr(u);
    }
    dsState_.sampleCount += 1;

    // Window = windowEddies eddy-passage times.
    const scalar tauEddy = max(dsState_.maxSigma, SMALL)/Ub;
    const scalar window = dsState_.windowEddies*tauEddy;
    if (t - dsState_.phaseStartTime >= window) finalizeDownstreamWindow();
}


void
Foam::turbulentSyntheticEddyInletFvPatchVectorField::finalizeDownstreamWindow()
{
    const scalar n = max(dsState_.sampleCount, scalar(1));
    const scalar Rscale = useUserR_ ? Rscale_ : sqr(Utau_)*Rscale_;

    symmTensor measWsum(symmTensor::zero);
    scalar areaSum = 0;
    forAll(dsState_.myProbes, pp)
    {
        const label p = dsState_.myProbes[pp];
        const vector Um = dsState_.sumU[p]/n;
        const symmTensor Rp = dsState_.sumUU[p]/n - sqr(Um);
        const tensor T = localFrameT(patchNormal_, dsState_.probeWn[p]);
        const symmTensor Rloc = globalToLocal(Rp, T)/max(Rscale, SMALL);
        const scalar a = dsState_.probeArea[p];
        measWsum += a*Rloc;
        areaSum += a;
    }
    reduce(measWsum, sumOp());
    reduce(areaSum, sumOp());
    areaSum = max(areaSum, SMALL);
    const symmTensor measMean = measWsum/areaSum;

    // Error from area-averaged TKE
    const symmTensor& rm = dsState_.refMean;
    const scalar trRef  = rm.xx() + rm.yy() + rm.zz();
    const scalar trMeas = measMean.xx() + measMean.yy() + measMean.zz();
    const scalar tkeErr = mag(trMeas - trRef)/max(trRef, SMALL);

    // tkeErr scores the current faceRlocal_. Track the best for the early stop.
    // A stall grows the next window x1.5 (cap 4x) to lower the noise floor.
    if (tkeErr < dsState_.bestTkeErr)
    {
        dsState_.bestTkeErr = tkeErr;
        dsState_.bestTarget = faceRlocal_;
        dsState_.noImprove = 0;
    }
    else
    {
        ++dsState_.noImprove;
        const scalar grown = min(1.5*dsState_.windowEddies, 4*ds_.windowEddies);
        if (grown > dsState_.windowEddies)
        {
            Info<< "  window grown " << dsState_.windowEddies << " -> " << grown
                << " eddy-passages (stall)" << endl;
            dsState_.windowEddies = grown;
        }
    }

    // Scale the whole target by one isotropic factor, under-relaxed (0.8) and clamped to [0.1, 10],
    // capped so no normal stress exceeds maxAmp*reference.
    const scalar relax = 0.8;
    scalar f = boundedScale(trRef, trMeas, relax, 0.1, 10.0);

    const direction nrm[3] = {0, 3, 5};   // xx, yy, zz
    for (label i = 0; i < 3; ++i)
    {
        const direction c = nrm[i];
        const scalar am = dsState_.appliedMean[c];
        if (am*f > ds_.maxAmp*rm[c] && am > SMALL) f = ds_.maxAmp*rm[c]/am;
    }

    forAll(faceRlocal_, i)
        if (faceActive_[i]) faceRlocal_[i] *= f;
    forAll(dsState_.probeTarget, p) dsState_.probeTarget[p] *= f;
    dsState_.appliedMean *= f;

    Info<< "  scaling f=" << f << " (uniform, all components)" << endl;
    applyDownstreamTarget();

    ++dsState_.pass;
    const scalar t = db().time().value();
    // Converged only once the TKE error stays below tol for 2 consecutive passes.
    if (tkeErr < ds_.tol) ++dsState_.belowTol; else dsState_.belowTol = 0;
    const bool converged = (dsState_.belowTol >= 2);
    const bool last = (dsState_.pass >= ds_.maxPasses);
    const bool stalled = (dsState_.noImprove >= ds_.patience);

    // TKE = 0.5*trace(R), in the controller's normalised form.
    Info<< "downstream R-correction pass " << dsState_.pass << " @t=" << t
        << ": window=" << label(dsState_.sampleCount) << " steps (~"
        << dsState_.windowEddies << " eddy-passages)  k=" << 0.5*trMeas
        << " (target " << 0.5*trRef << ")  TKE err=" << tkeErr
        << "  (best=" << dsState_.bestTkeErr << ")"
        << (converged ? "  [converged]"
              : (stalled ? "  [stalled]" : (last ? "  [max passes]" : "")))
        << endl;

    if (converged || last || stalled)
    {
        dsState_.ph = downstreamState::phase::done;
        // End on the best-TKE-error target seen
        if (dsState_.bestTarget.size() == faceRlocal_.size())
        {
            faceRlocal_ = dsState_.bestTarget;
            applyDownstreamTarget();
            Info<< "  reverted to best target (TKE err=" << dsState_.bestTkeErr
                << ")" << endl;
        }
        writeConvergedR();
    }
    else
    {
        dsState_.ph = downstreamState::phase::warmup;
        dsState_.phaseStartTime = t;
    }
}


void
Foam::turbulentSyntheticEddyInletFvPatchVectorField::applyDownstreamTarget()
{
    rebuildTemplates();
    gatherFaceData();
}


void Foam::turbulentSyntheticEddyInletFvPatchVectorField::writeConvergedR() const
{
    // Write the converged inlet R to the startTime directory as "R_converged",
    // in the global frame and physical units (the form a user R is read in). A
    // similar case can restart from it via customRName=R_converged.
    const fvMesh& mesh = this->internalField().mesh();
    const word t0 = mesh.time().timeName(mesh.time().startTime().value());
    // register=false: transient, write-only field.
    volSymmTensorField Rc
    (
        IOobject
        (
            "R_converged", t0, mesh,
            IOobject::NO_READ, IOobject::NO_WRITE, false
        ),
        mesh,
        dimensionedSymmTensor(sqr(dimVelocity), Zero)
    );
    const label pIdx = patch().index();
    const scalar physScale = sqr(Utau_)*Rscale_;
    forAll(faceRlocal_, i)
        Rc.boundaryFieldRef()[pIdx][i] =
            localToGlobal(faceRlocal_[i], localFrameT(patchNormal_, wallNormal_[i]))
          * physScale;
    // Force the startTime directory.
    Rc.instance() = t0;
    fileHandler().writeObject
    (
        Rc, mesh.time().writeFormat(), IOstream::currentVersion,
        mesh.time().writeCompression()
    );
    Info<< "  wrote converged inlet R -> " << t0 << "/R_converged" << endl;
}


// * * * * * * * * * Ghost-transform construction * * * * * * * * * * * * * //

void Foam::turbulentSyntheticEddyInletFvPatchVectorField::buildGhostTransforms()
{
    cyclicTranslations_.clear();

    const polyBoundaryMesh& bm =
        patch().poly().boundaryMesh();

    DynamicList<word> seen;
    forAll(bm, patchI)
    {
        if (!isA<cyclicPolyPatch>(bm[patchI]))
        {
            continue;
        }
        const cyclicPolyPatch& cpp =
            refCast<const cyclicPolyPatch>(bm[patchI]);
        bool already = false;
        forAll(seen, s)
        {
            if (seen[s] == cpp.name()) { already = true; break; }
        }
        if (already)
        {
            continue;
        }
        const cyclicPolyPatch& nbr = cpp.nbrPatch();
        seen.append(cpp.name());
        seen.append(nbr.name());

        // Pull translation from the cyclic-patch metadata.
        vector t = Zero;
        if (cpp.transformComplete())
        {
            t = cpp.transform().t();
        }
        else
        {
            // Fallback: average face centres if the metadata is not
            // populated. Skip if globally empty.
            const vectorField cppCf(cpp.faceCentres());
            const vectorField nbrCf(nbr.faceCentres());
            const label nGlobal =
                returnReduce(cppCf.size(), sumOp());
            if (nGlobal == 0)
            {
                continue;
            }
            t = gAverage(nbrCf) - gAverage(cppCf);
        }
        if (mag(t) > SMALL)
        {
            cyclicTranslations_.append(t);
        }
    }
    if (debug)
    {
        Info<< "    BC '" << patch().name() << "' built "
            << cyclicTranslations_.size() << " cyclic translations:";
        forAll(cyclicTranslations_, c)
        {
            Info<< " t" << c << "=" << cyclicTranslations_[c];
        }
        Info<< endl;
    }
}


// * * * * * * * * * * * * * Parallel gather of face data * * * * * * * * * //

void Foam::turbulentSyntheticEddyInletFvPatchVectorField::gatherFaceData()
{
    // Gather every face onto every rank (concatenated in proc order).
    gCf_ = gatherGlobal<vector>(patch().Cf());
    gMagSf_ = gatherGlobal<scalar>(patch().magSf());
    gSigmaX_ = gatherGlobal<scalar>(faceSigmaX_);
    gRpg_ = gatherGlobal<tensor>(faceRpg_);
    gGamma_ = gatherGlobal<scalar>(faceGamma_);
    gAlphaCoeff_ = gatherGlobal<vector>(faceAlphaCoeff_);
    gActive_ = gatherGlobal<bool>(faceActive_);

    // Area-weighted spawn CDF (eddy density proportional to face area).
    const label N = gMagSf_.size();
    gCDF_.setSize(N);
    scalar cum = 0;
    forAll(gMagSf_, i)
    {
        cum += gMagSf_[i];
        gCDF_[i] = cum;
    }
}


// * * * * * * * * * * * * * * * Eddy spawning * * * * * * * * * * * * * * * //

Foam::syntheticEddy
Foam::turbulentSyntheticEddyInletFvPatchVectorField::oneEddy(scalar depth)
{
    if (gCDF_.empty() || gCDF_.last() <= SMALL)
    {
        FatalErrorInFunction
            << "Patch '" << patch().name() << "' has no faces with "
            << "positive area in the global CDF (gCDF size="
            << gCDF_.size()
            << ", total weight=" << (gCDF_.empty() ? 0 : gCDF_.last())
            << "); cannot spawn eddies."
            << exit(FatalError);
    }

    const scalar u = randGen_.scalar01()*gCDF_.last();
    label j = upperBound(gCDF_.begin(), gCDF_.size(), u);
    if (j >= gCDF_.size())
    {
        j = gCDF_.size() - 1;
    }

    const scalar e0 = randGen_.scalar01() < 0.5 ? 1.0 : -1.0;
    const scalar e1 = randGen_.scalar01() < 0.5 ? 1.0 : -1.0;
    const scalar e2 = randGen_.scalar01() < 0.5 ? 1.0 : -1.0;

    // Per-eddy convection velocity.
    const scalar Uconv = Umean_;

    if (!gActive_[j])
    {
        return syntheticEddy(gCf_[j], depth, vector::zero, vector::zero, tensor::I, -1, Uconv);
    }

    const scalar sigmaX = gSigmaX_[j];
    if (sigmaX <= 0)
    {
        return syntheticEddy(gCf_[j], depth, vector::zero, vector::zero, tensor::I, -1, Uconv);
    }

    const scalar g = gGamma_[j];
    const vector sigma(sigmaX/g, sigmaX/g, sigmaX);

    // Dimensionless shape amplitude (scaled to physical by c2 in updateCoeffs).
    const scalar coeff = 1.0/sigmaX;
    const vector alpha
    (
        e0*gAlphaCoeff_[j][0]*coeff,
        e1*gAlphaCoeff_[j][1]*coeff,
        e2*gAlphaCoeff_[j][2]*coeff
    );

    const scalar c1 =
        g/sqrt(Foam::constant::mathematical::pi*4.0/3.0*pow3(sigmaX));

    return syntheticEddy(gCf_[j], depth, sigma, alpha, gRpg_[j], c1, Uconv);
}


void Foam::turbulentSyntheticEddyInletFvPatchVectorField::initializeEddies()
{
    DynamicList<syntheticEddy> local(size());
    scalar sumVol = 0;
    while (sumVol/v0_ < eddyDensity_)
    {
        const scalar depth = randGen_.scalarAB(-maxSigmaX_, maxSigmaX_);
        syntheticEddy e = oneEddy(depth);
        local.append(e);
        sumVol += e.volume();
    }
    eddies_.transfer(local);
    nEddy_ = eddies_.size();

    if (debug)
        Info<< "Initialised patch " << patch().name() << " with "
            << nEddy_ << " eddies (v0=" << v0_
            << ", maxSigmaX=" << maxSigmaX_ << ")." << endl;
}


void Foam::turbulentSyntheticEddyInletFvPatchVectorField::convectEddies()
{
    const scalar deltaT = db().time().deltaTValue();

    updateUmean();
    forAll(eddies_, eI)
    {
        syntheticEddy& e = eddies_[eI];
        e.move(deltaT*e.Uconv()*eddyConvectionFactor_);
        const scalar x = e.x();
        if (x > maxSigmaX_)
        {
            e = oneEddy(x - 2*floor(x/maxSigmaX_)*maxSigmaX_);
        }
    }
}


// * * * * * * * * * * * uPrime evaluation at a point * * * * * * * * * * * //

Foam::vector
Foam::turbulentSyntheticEddyInletFvPatchVectorField::uPrimeEddy(const vector& pos) const
{
    vector uPrime(Zero);

    forAll(eddies_, k)
    {
        const syntheticEddy& e = eddies_[k];

        uPrime += e.uPrime(pos, patchNormal_);

        forAll(cyclicTranslations_, c)
        {
            const vector& t = cyclicTranslations_[c];
            vector shifted = pos + t;
            if (mag(shifted - e.position0()) < maxSigmaX_)
            {
                uPrime += e.uPrime(shifted, patchNormal_);
            }
            shifted = pos - t;
            if (mag(shifted - e.position0()) < maxSigmaX_)
            {
                uPrime += e.uPrime(shifted, patchNormal_);
            }
        }
    }
    return uPrime;
}

// * * * * * * * * * * * * * * * updateCoeffs * * * * * * * * * * * * * * * //

void Foam::turbulentSyntheticEddyInletFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    if (!initialized_)
    {
        initializeGeometry();
        buildGhostTransforms();
        initializeProfiles();
        gatherFaceData();
        initializeEddies();
        initialized_ = true;
    }

    convectEddies();

    updateUmean();

    // Quasi-steady Utau tracking: re-evaluate from the live Umean
    if
    (
        !UtauIsInDict_
     && (ds_.distance <= 0 || dsState_.ph != downstreamState::phase::accumulate)
    )
    {
        const scalar UtauNew = frictionUtau(Umean_, avNu_);
        if (mag(UtauNew/max(Utau_, SMALL) - 1) > 1e-4)
        {
            if (debug)
                Info<< "Utau tracking @t=" << db().time().value()
                    << ": Umean=" << Umean_ << "  Utau " << Utau_
                    << " -> " << UtauNew << endl;
            Utau_ = UtauNew;
        }
    }

    // Opt-in downstream R-correction
    stepDownstreamController();

    const scalar c2 = sqrt((315.0/16.0)*v0_*Rscale_/scalar(nEddy_))
                    * (useUserR_ ? 1.0 : Utau_);

    const pointField& Cf = patch().Cf();
    vectorField uEddy(patch().size(), Zero);
    forAll(uEddy, faceI)
    {
        uEddy[faceI] = c2*uPrimeEddy(Cf[faceI]);
    }

    {
        const scalarField magSf(patch().magSf());
        const scalar Atot = max(gSum(magSf), SMALL);
        const vector meanEddy
        (
            gSum(magSf*uEddy.component(0))/Atot,
            gSum(magSf*uEddy.component(1))/Atot,
            gSum(magSf*uEddy.component(2))/Atot
        );
        uEddy -= meanEddy;
    }

    vectorField Upatch = U_*Umean_ + uEddy;

    vectorField::operator=(Upatch);
    fixedValueFvPatchVectorField::updateCoeffs();
}


// * * * * * * * * * * * * * * * write * * * * * * * * * * * * * * * * * * //

void Foam::turbulentSyntheticEddyInletFvPatchVectorField::write(Ostream& os) const
{
    fvPatchVectorField::write(os);
    writeEntry
    (
        os,
        this->db().time().userUnits(),
        this->internalField().dimensions(),
        inFlow_()
    );
    writeEntryIfDifferent<scalar>(os, "eddyConvectionFactor", 0.5, eddyConvectionFactor_);
    if (UtauIsInDict_) writeEntry(os, "Utau", Utau_);
    writeEntryIfDifferent<scalar>(os, "swirlNumber", 0, swirlNumber_);
    if (swirlAxisIsInDict_)   writeEntry(os, "swirlAxis", swirlAxis_);
    if (swirlCenterIsInDict_) writeEntry(os, "swirlCenter", swirlCenter_);
    if (wallPatchesIsInDict_) writeEntry(os, "wallPatches", wallPatchNames_);
    writeEntryIfDifferent<bool>(os, "inletRcorrection", false, inletRcorrection_);
    if (inletRcorrection_)
        writeEntryIfDifferent<scalar>
            (os, "inletRcorrectionOverlaps", 500, mcOverlaps_);
    if (!customRName_.empty())     writeEntry(os, "customRName", customRName_);
    if (!customLName_.empty())     writeEntry(os, "customLName", customLName_);
    if (!customUmeanName_.empty()) writeEntry(os, "customUmeanName", customUmeanName_);
    writeEntryIfDifferent<scalar>(os, "eddyDensity", 5, eddyDensity_);
    writeEntryIfDifferent<scalar>(os, "Lscale", 1, Lscale_);
    writeEntryIfDifferent<scalar>(os, "Rscale", 1, Rscale_);
    writeEntryIfDifferent<label>(os, "nCellPerEddy", 3, nCellPerEddy_);
    writeEntryIfDifferent<scalar>(os, "rhoInf", 1, rhoInf_);
    writeEntryIfDifferent<label>(os, "seed", 536, seedInput_);

    // Emit downstream settings only when enabled.
    if (ds_.distance > 0)
    {
        writeEntry(os, "correctionDistance", ds_.distance);
        writeEntryIfDifferent<scalar>
            (os, "downstreamWindowEddies", 10, ds_.windowEddies);
        writeEntryIfDifferent<label>
            (os, "downstreamMaxPasses", 10, ds_.maxPasses);
        writeEntryIfDifferent<scalar>
            (os, "downstreamStartTime", 0, ds_.startTime);
        writeEntryIfDifferent<scalar>
            (os, "downstreamMaxAmp", 10, ds_.maxAmp);
        writeEntryIfDifferent<label>
            (os, "downstreamPatience", 3, ds_.patience);
        writeEntryIfDifferent<scalar>
            (os, "downstreamTol", 0.05, ds_.tol);
    }
}


// * * * * * * * * * * * * * * * Run-time selection * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchVectorField,
        turbulentSyntheticEddyInletFvPatchVectorField
    );
}


// ************************************************************************* //
