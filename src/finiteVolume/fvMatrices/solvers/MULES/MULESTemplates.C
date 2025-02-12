/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2025 OpenFOAM Foundation
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

#include "MULES.H"
#include "upwind.H"
#include "fvcSurfaceIntegrate.H"
#include "localEulerDdtScheme.H"
#include "wedgeFvPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class RdeltaTType, class RhoType, class SpType, class SuType>
void Foam::MULES::explicitSolve
(
    const RdeltaTType& rDeltaT,
    const RhoType& rho,
    volScalarField& psi,
    const surfaceScalarField& phiPsi,
    const SpType& Sp,
    const SuType& Su
)
{
    Info<< "MULES: Solving for " << psi.name() << endl;

    const fvMesh& mesh = psi.mesh();

    scalarField& psiIf = psi;
    const scalarField& psi0 = psi.oldTime();

    psiIf = 0.0;
    fvc::surfaceIntegrate(psiIf, phiPsi);

    if (mesh.moving())
    {
        psiIf =
        (
            mesh.Vsc0()().primitiveField()*rho.oldTime().primitiveField()
           *psi0*rDeltaT/mesh.Vsc()().primitiveField()
          + Su.primitiveField()
          - psiIf
        )/(rho.primitiveField()*rDeltaT - Sp.primitiveField());
    }
    else
    {
        psiIf =
        (
            rho.oldTime().primitiveField()*psi0*rDeltaT
          + Su.primitiveField()
          - psiIf
        )/(rho.primitiveField()*rDeltaT - Sp.primitiveField());
    }

    psi.correctBoundaryConditions();
}


template<class RhoType>
void Foam::MULES::explicitSolve
(
    const RhoType& rho,
    volScalarField& psi,
    const surfaceScalarField& phiPsi
)
{
    explicitSolve(rho, psi, phiPsi, zeroField(), zeroField());
}


template<class RhoType, class SpType, class SuType>
void Foam::MULES::explicitSolve
(
    const RhoType& rho,
    volScalarField& psi,
    const surfaceScalarField& phiPsi,
    const SpType& Sp,
    const SuType& Su
)
{
    const fvMesh& mesh = psi.mesh();

    if (fv::localEulerDdt::enabled(mesh))
    {
        const volScalarField& rDeltaT = fv::localEulerDdt::localRDeltaT(mesh);
        explicitSolve(rDeltaT, rho, psi, phiPsi, Sp, Su);
    }
    else
    {
        const scalar rDeltaT = 1.0/mesh.time().deltaTValue();
        explicitSolve(rDeltaT, rho, psi, phiPsi, Sp, Su);
    }
}


template<class RhoType, class PsiMaxType, class PsiMinType>
void Foam::MULES::explicitSolve
(
    const RhoType& rho,
    volScalarField& psi,
    const surfaceScalarField& phiBD,
    surfaceScalarField& phiPsi,
    const PsiMaxType& psiMax,
    const PsiMinType& psiMin
)
{
    explicitSolve
    (
        rho,
        psi,
        phiBD,
        phiPsi,
        zeroField(),
        zeroField(),
        psiMax,
        psiMin
    );
}


template
<
    class RhoType,
    class SpType,
    class SuType,
    class PsiMaxType,
    class PsiMinType
>
void Foam::MULES::explicitSolve
(
    const RhoType& rho,
    volScalarField& psi,
    const surfaceScalarField& phi,
    surfaceScalarField& phiPsi,
    const SpType& Sp,
    const SuType& Su,
    const PsiMaxType& psiMax,
    const PsiMinType& psiMin
)
{
    const fvMesh& mesh = psi.mesh();

    psi.correctBoundaryConditions();

    if (fv::localEulerDdt::enabled(mesh))
    {
        const volScalarField& rDeltaT = fv::localEulerDdt::localRDeltaT(mesh);
        limit(rDeltaT, rho, psi, phi, phiPsi, Sp, Su, psiMax, psiMin, false);
        explicitSolve(rDeltaT, rho, psi, phiPsi, Sp, Su);
    }
    else
    {
        const scalar rDeltaT = 1.0/mesh.time().deltaTValue();
        limit(rDeltaT, rho, psi, phi, phiPsi, Sp, Su, psiMax, psiMin, false);
        explicitSolve(rDeltaT, rho, psi, phiPsi, Sp, Su);
    }
}


template
<
    class RdeltaTType,
    class RhoType,
    class SpType,
    class PsiMaxType,
    class PsiMinType
>
void Foam::MULES::limiter
(
    surfaceScalarField& lambda,
    const RdeltaTType& rDeltaT,
    const RhoType& rho,
    const volScalarField& psi,
    const scalarField& SuCorr,
    const surfaceScalarField& phiBD,
    const surfaceScalarField& phiCorr,
    const SpType& Sp,
    const PsiMaxType& psiMax,
    const PsiMinType& psiMin
)
{
    const scalarField& psiIf = psi;
    const volScalarField::Boundary& psiBf = psi.boundaryField();

    const fvMesh& mesh = psi.mesh();

    const dictionary& MULEScontrols = mesh.solution().solverDict(psi.name());

    const label nLimiterIter
    (
        MULEScontrols.lookupOrDefault<label>("nLimiterIter", 3)
    );

    const scalar smoothLimiter
    (
        MULEScontrols.lookupOrDefault<scalar>("smoothLimiter", 0)
    );

    const scalar extremaCoeff
    (
        MULEScontrols.lookupOrDefault<scalar>("extremaCoeff", 0)
    );

    const scalar boundaryExtremaCoeff
    (
        MULEScontrols.lookupOrDefault<scalar>
        (
            "boundaryExtremaCoeff",
            extremaCoeff
        )
    );

    const scalar boundaryDeltaExtremaCoeff
    (
        max(boundaryExtremaCoeff - extremaCoeff, 0)
    );

    const scalar tol
    (
        nLimiterIter == 1
      ? 0
      : MULEScontrols.lookupOrDefault<scalar>("MULEStolerance", 0)
    );

    const labelUList& owner = mesh.owner();
    const labelUList& neighb = mesh.neighbour();
    tmp<volScalarField::Internal> tVsc = mesh.Vsc();
    const scalarField& V = tVsc();

    const surfaceScalarField::Boundary& phiBDBf = phiBD.boundaryField();

    const scalarField& phiCorrIf = phiCorr;
    const surfaceScalarField::Boundary& phiCorrBf = phiCorr.boundaryField();

    scalarField& lambdaIf = lambda;
    surfaceScalarField::Boundary& lambdaBf = lambda.boundaryFieldRef();

    scalarField psiMaxn(psiIf.size());
    scalarField psiMinn(psiIf.size());

    psiMaxn = psiMin;
    psiMinn = psiMax;

    scalarField phiCorrNorm;
    if (tol != 0)
    {
        phiCorrNorm = (V*(rho.primitiveField()*rDeltaT - Sp.primitiveField()));
    }

    scalarField sumPhip(psiIf.size(), 0.0);
    scalarField mSumPhim(psiIf.size(), 0.0);

    forAll(phiCorrIf, facei)
    {
        const label own = owner[facei];
        const label nei = neighb[facei];

        psiMaxn[own] = max(psiMaxn[own], psiIf[nei]);
        psiMinn[own] = min(psiMinn[own], psiIf[nei]);

        psiMaxn[nei] = max(psiMaxn[nei], psiIf[own]);
        psiMinn[nei] = min(psiMinn[nei], psiIf[own]);

        const scalar phiCorrf = phiCorrIf[facei];

        if (phiCorrf > 0)
        {
            sumPhip[own] += phiCorrf;
            mSumPhim[nei] += phiCorrf;
        }
        else
        {
            mSumPhim[own] -= phiCorrf;
            sumPhip[nei] -= phiCorrf;
        }
    }

    forAll(phiCorrBf, patchi)
    {
        const fvPatchScalarField& psiPf = psiBf[patchi];
        const scalarField& phiCorrPf = phiCorrBf[patchi];

        const labelList& pFaceCells = mesh.boundary()[patchi].faceCells();

        if (psiPf.coupled())
        {
            const scalarField psiPNf(psiPf.patchNeighbourField());

            forAll(phiCorrPf, pFacei)
            {
                const label pfCelli = pFaceCells[pFacei];

                psiMaxn[pfCelli] = max(psiMaxn[pfCelli], psiPNf[pFacei]);
                psiMinn[pfCelli] = min(psiMinn[pfCelli], psiPNf[pFacei]);
            }
        }
        else if (psiPf.fixesValue())
        {
            forAll(phiCorrPf, pFacei)
            {
                const label pfCelli = pFaceCells[pFacei];

                psiMaxn[pfCelli] = max(psiMaxn[pfCelli], psiPf[pFacei]);
                psiMinn[pfCelli] = min(psiMinn[pfCelli], psiPf[pFacei]);
            }
        }
        else
        {
            // Add the optional additional allowed boundary extrema
            if (boundaryDeltaExtremaCoeff > 0)
            {
                forAll(phiCorrPf, pFacei)
                {
                    const label pfCelli = pFaceCells[pFacei];

                    const scalar extrema =
                        boundaryDeltaExtremaCoeff
                       *(psiMax[pfCelli] - psiMin[pfCelli]);

                    psiMaxn[pfCelli] += extrema;
                    psiMinn[pfCelli] -= extrema;
                }
            }
        }

        forAll(phiCorrPf, pFacei)
        {
            const label pfCelli = pFaceCells[pFacei];
            const scalar phiCorrf = phiCorrPf[pFacei];

            if (phiCorrf > 0)
            {
                sumPhip[pfCelli] += phiCorrf;
            }
            else
            {
                mSumPhim[pfCelli] -= phiCorrf;
            }
        }
    }

    if (extremaCoeff > 0)
    {
        psiMaxn = min(psiMaxn + extremaCoeff*(psiMax - psiMin), psiMax);
        psiMinn = max(psiMinn - extremaCoeff*(psiMax - psiMin), psiMin);
    }
    else
    {
        psiMaxn = min(psiMaxn, psiMax);
        psiMinn = max(psiMinn, psiMin);
    }

    if (smoothLimiter > small)
    {
        psiMaxn =
            min(smoothLimiter*psiIf + (1.0 - smoothLimiter)*psiMaxn, psiMax);
        psiMinn =
            max(smoothLimiter*psiIf + (1.0 - smoothLimiter)*psiMinn, psiMin);
    }

    psiMaxn =
        V*((rho.primitiveField()*rDeltaT - Sp.primitiveField())*psiMaxn)
      - SuCorr;

    psiMinn =
        SuCorr
      - V*((rho.primitiveField()*rDeltaT - Sp.primitiveField())*psiMinn);

    scalarField sumlPhip(psiIf.size());
    scalarField mSumlPhim(psiIf.size());

    // Allocate storage for lambda0 on coupled patches
    // for optional convergence test
    surfaceScalarField::Boundary lambdaBf0(mesh.boundary());
    if (tol != 0)
    {
        forAll(lambdaBf, patchi)
        {
            fvsPatchScalarField& lambdaPf = lambdaBf[patchi];

            if (lambdaPf.coupled())
            {
                lambdaBf0.set
                (
                    patchi,
                    new calculatedFvsPatchField<scalar>
                    (
                        mesh.boundary()[patchi],
                        surfaceScalarField::Internal::null()
                    )
                );
            }
        }
    }

    for (int j=0; j<nLimiterIter; j++)
    {
        // Convergence test parameter
        scalar maxDeltaLambdaPhiCorrRes = 0;

        // Sum limited positive and negative fluxes
        // Not needed for first iteration
        if (j > 0)
        {
            sumlPhip = 0;
            mSumlPhim = 0;

            forAll(lambdaIf, facei)
            {
                const label own = owner[facei];
                const label nei = neighb[facei];

                const scalar lambdaPhiCorrf = lambdaIf[facei]*phiCorrIf[facei];

                if (lambdaPhiCorrf > 0)
                {
                    sumlPhip[own] += lambdaPhiCorrf;
                    mSumlPhim[nei] += lambdaPhiCorrf;
                }
                else
                {
                    mSumlPhim[own] -= lambdaPhiCorrf;
                    sumlPhip[nei] -= lambdaPhiCorrf;
                }
            }

            forAll(lambdaBf, patchi)
            {
                scalarField& lambdaPf = lambdaBf[patchi];
                const scalarField& phiCorrfPf = phiCorrBf[patchi];

                const labelList& pFaceCells =
                    mesh.boundary()[patchi].faceCells();

                forAll(lambdaPf, pFacei)
                {
                    const label pfCelli = pFaceCells[pFacei];
                    const scalar lambdaPhiCorrf =
                        lambdaPf[pFacei]*phiCorrfPf[pFacei];

                    if (lambdaPhiCorrf > 0)
                    {
                        sumlPhip[pfCelli] += lambdaPhiCorrf;
                    }
                    else
                    {
                        mSumlPhim[pfCelli] -= lambdaPhiCorrf;
                    }
                }
            }
        }

        // Reuse storage of sumlPhip and mSumlPhim for lambdam and lambdap
        scalarField& lambdam = sumlPhip;
        scalarField& lambdap = mSumlPhim;

        if (j == 0)
        {
            forAll(lambdam, celli)
            {
                lambdam[celli] =
                    max(min
                    (
                        psiMaxn[celli]/(mSumPhim[celli] + rootVSmall),
                        1.0), 0.0
                    );

                lambdap[celli] =
                    max(min
                    (
                        psiMinn[celli]/(sumPhip[celli] + rootVSmall),
                        1.0), 0.0
                    );
            }
        }
        else
        {
            forAll(lambdam, celli)
            {
                lambdam[celli] =
                    max(min
                    (
                        (sumlPhip[celli] + psiMaxn[celli])
                       /(mSumPhim[celli] + rootVSmall),
                        1.0), 0.0
                    );

                lambdap[celli] =
                    max(min
                    (
                        (mSumlPhim[celli] + psiMinn[celli])
                       /(sumPhip[celli] + rootVSmall),
                        1.0), 0.0
                    );
            }
        }

        forAll(lambdaIf, facei)
        {
            const scalar lambdaIf0 = lambdaIf[facei];

            if (phiCorrIf[facei] > 0)
            {
                lambdaIf[facei] =
                    min(lambdap[owner[facei]], lambdam[neighb[facei]]);
            }
            else
            {
                lambdaIf[facei] =
                    min(lambdam[owner[facei]], lambdap[neighb[facei]]);
            }

            if (tol > 0)
            {
                const scalar phiCorrRes =
                    mag(phiCorrIf[facei])
                   /min(phiCorrNorm[owner[facei]], phiCorrNorm[neighb[facei]]);

                if (phiCorrRes > tol)
                {
                    maxDeltaLambdaPhiCorrRes = max
                    (
                        maxDeltaLambdaPhiCorrRes,
                        mag(lambdaIf[facei] - lambdaIf0)*phiCorrRes
                    );
                }
            }
        }

        // Take minimum of value across coupled patches
        forAll(lambdaBf, patchi)
        {
            fvsPatchScalarField& lambdaPf = lambdaBf[patchi];
            const scalarField& phiCorrfPf = phiCorrBf[patchi];
            const fvPatchScalarField& psiPf = psiBf[patchi];

            if (isA<wedgeFvPatch>(mesh.boundary()[patchi]))
            {
                lambdaPf = 0;
            }
            else if (psiPf.coupled())
            {
                const labelList& pFaceCells =
                    mesh.boundary()[patchi].faceCells();

                if (tol > 0)
                {
                    lambdaBf0[patchi] = lambdaPf;
                }

                forAll(lambdaPf, pFacei)
                {
                    const label pfCelli = pFaceCells[pFacei];

                    if (phiCorrfPf[pFacei] > 0)
                    {
                        lambdaPf[pFacei] = lambdap[pfCelli];
                    }
                    else
                    {
                        lambdaPf[pFacei] = lambdam[pfCelli];
                    }
                }
            }
            else
            {
                const labelList& pFaceCells =
                    mesh.boundary()[patchi].faceCells();

                const scalarField& phiBDPf = phiBDBf[patchi];

                forAll(lambdaPf, pFacei)
                {
                    // Limit outlet faces only
                    if ((phiBDPf[pFacei] + phiCorrfPf[pFacei]) > small*small)
                    {
                        const label pfCelli = pFaceCells[pFacei];
                        const scalar lambdaPf0 = lambdaPf[pFacei];

                        if (phiCorrfPf[pFacei] > 0)
                        {
                            lambdaPf[pFacei] = lambdap[pfCelli];
                        }
                        else
                        {
                            lambdaPf[pFacei] = lambdam[pfCelli];
                        }

                        if (tol > 0)
                        {
                            const scalar phiCorrRes =
                                mag(phiCorrfPf[pFacei])/phiCorrNorm[pfCelli];

                            if (phiCorrRes > tol)
                            {
                                maxDeltaLambdaPhiCorrRes = max
                                (
                                    maxDeltaLambdaPhiCorrRes,
                                    mag(lambdaPf[pFacei] - lambdaPf0)*phiCorrRes
                                );
                            }
                        }
                    }
                }
            }
        }

        // Take minimum value of limiter across coupled patches
        surfaceScalarField::Boundary lambdaNbrBf
        (
            surfaceScalarField::Internal::null(),
            lambdaBf.boundaryNeighbourField()
        );

        forAll(lambdaBf, patchi)
        {
            fvsPatchScalarField& lambdaPf = lambdaBf[patchi];

            if (lambdaPf.coupled())
            {
                const fvsPatchScalarField& lambdaNbrPf = lambdaNbrBf[patchi];
                lambdaPf = min(lambdaPf, lambdaNbrPf);

                if (tol > 0)
                {
                    const fvsPatchScalarField& lambdaPf0 = lambdaBf0[patchi];
                    const scalarField& phiCorrfPf = phiCorrBf[patchi];

                    const labelList& pFaceCells =
                        mesh.boundary()[patchi].faceCells();

                    forAll(lambdaPf, pFacei)
                    {
                        const scalar phiCorrRes =
                            mag(phiCorrfPf[pFacei])
                           /phiCorrNorm[pFaceCells[pFacei]];

                        if (phiCorrRes > tol)
                        {
                            maxDeltaLambdaPhiCorrRes = max
                            (
                                maxDeltaLambdaPhiCorrRes,
                                mag(lambdaPf[pFacei] - lambdaPf0[pFacei])
                               *phiCorrRes
                            );
                        }
                    }
                }
            }
        }

        // Optional convergence test
        if (tol != 0)
        {
            reduce(maxDeltaLambdaPhiCorrRes, maxOp<scalar>());

            if (debug)
            {
                Info<< "MULES: maxDeltaLambdaPhiCorrRes "
                    << maxDeltaLambdaPhiCorrRes << endl;
            }

            if (maxDeltaLambdaPhiCorrRes < tol) break;
        }
    }
}


template
<
    class RdeltaTType,
    class RhoType,
    class SpType,
    class SuType,
    class PsiMaxType,
    class PsiMinType
>
void Foam::MULES::limit
(
    const RdeltaTType& rDeltaT,
    const RhoType& rho,
    const volScalarField& psi,
    const surfaceScalarField& phi,
    surfaceScalarField& phiPsi,
    const SpType& Sp,
    const SuType& Su,
    const PsiMaxType& psiMax,
    const PsiMinType& psiMin,
    const bool returnCorr
)
{
    const fvMesh& mesh = psi.mesh();

    surfaceScalarField phiBD(upwind<scalar>(psi.mesh(), phi).flux(psi));

    const scalarField& phiBDIf = phiBD;
    surfaceScalarField::Boundary& phiBDBf = phiBD.boundaryFieldRef();

    const surfaceScalarField::Boundary& phiPsiBf = phiPsi.boundaryField();

    const labelUList& owner = mesh.owner();
    const labelUList& neighb = mesh.neighbour();

    forAll(phiBDBf, patchi)
    {
        fvsPatchScalarField& phiBDPf = phiBDBf[patchi];

        if (!phiBDPf.coupled())
        {
            phiBDPf = phiPsiBf[patchi];
        }
    }

    surfaceScalarField& phiCorr = phiPsi;
    phiCorr -= phiBD;

    tmp<volScalarField::Internal> tVsc = mesh.Vsc();
    const scalarField& V = tVsc();

    // Correction equation source
    scalarField SuCorr
    (
        mesh.moving()
      ? (mesh.Vsc0()().primitiveField()*rDeltaT*rho.oldTime().primitiveField())
       *psi.oldTime().primitiveField()
      + V*Su.primitiveField()
      : V
       *(
            (rho.oldTime().primitiveField()*rDeltaT)
           *psi.oldTime().primitiveField()
          + Su.primitiveField()
        )
    );

    // Subtract the sum of the bounded fluxes
    // from the correction equation source
    forAll(phiBDIf, facei)
    {
        SuCorr[owner[facei]] -= phiBDIf[facei];
        SuCorr[neighb[facei]] += phiBDIf[facei];
    }

    forAll(phiBDBf, patchi)
    {
        const scalarField& phiBDPf = phiBDBf[patchi];
        const labelList& pFaceCells = mesh.boundary()[patchi].faceCells();

        forAll(phiBDPf, pFacei)
        {
            SuCorr[pFaceCells[pFacei]] -= phiBDPf[pFacei];
        }
    }

    surfaceScalarField lambda
    (
        IOobject
        (
            "lambda",
            mesh.time().name(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        mesh,
        dimensionedScalar(dimless, 1)
    );

    limiter
    (
        lambda,
        rDeltaT,
        rho,
        psi,
        SuCorr,
        phiBD,
        phiCorr,
        Sp,
        psiMax,
        psiMin
    );

    if (returnCorr)
    {
        phiCorr *= lambda;
    }
    else
    {
        phiPsi = phiBD + lambda*phiCorr;
    }
}


template
<
    class RhoType,
    class SpType,
    class SuType,
    class PsiMaxType,
    class PsiMinType
>
void Foam::MULES::limit
(
    const RhoType& rho,
    const volScalarField& psi,
    const surfaceScalarField& phi,
    surfaceScalarField& phiPsi,
    const SpType& Sp,
    const SuType& Su,
    const PsiMaxType& psiMax,
    const PsiMinType& psiMin,
    const bool rtnCorr
)
{
    const fvMesh& mesh = psi.mesh();

    if (fv::localEulerDdt::enabled(mesh))
    {
        const volScalarField& rDeltaT = fv::localEulerDdt::localRDeltaT(mesh);
        limit(rDeltaT, rho, psi, phi, phiPsi, Sp, Su, psiMax, psiMin, rtnCorr);
    }
    else
    {
        const scalar rDeltaT = 1.0/mesh.time().deltaTValue();
        limit(rDeltaT, rho, psi, phi, phiPsi, Sp, Su, psiMax, psiMin, rtnCorr);
    }
}


template<template<class> class AlphaList, template<class> class PhiList>
void Foam::MULES::limitSum
(
    const AlphaList<const volScalarField>& alphas,
    PhiList<surfaceScalarField>& phiPsis,
    const surfaceScalarField& phi
)
{
    PtrList<surfaceScalarField> alphaPhiUDs(phiPsis.size());

    forAll(phiPsis, phasei)
    {
        alphaPhiUDs.set
        (
            phasei,
            upwind<scalar>(phi.mesh(), phi).flux(alphas[phasei])
        );

        phiPsis[phasei] -= alphaPhiUDs[phasei];
    }

    {
        UPtrList<scalarField> phiPsisInternal(phiPsis.size());

        forAll(phiPsisInternal, phasei)
        {
            phiPsisInternal.set(phasei, &phiPsis[phasei]);
        }

        limitSum(phiPsisInternal);
    }

    const surfaceScalarField::Boundary& phibf = phi.boundaryField();

    forAll(phibf, patchi)
    {
        if (phibf[patchi].coupled())
        {
            UPtrList<scalarField> phiPsisPatch(phiPsis.size());

            forAll(phiPsisPatch, phasei)
            {
                phiPsisPatch.set
                (
                    phasei,
                    &phiPsis[phasei].boundaryFieldRef()[patchi]
                );
            }

            limitSum(phiPsisPatch);
        }
    }

    forAll(phiPsis, phasei)
    {
        phiPsis[phasei] += alphaPhiUDs[phasei];
    }
}


// ************************************************************************* //
