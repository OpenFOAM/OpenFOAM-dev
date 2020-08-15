/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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
#include "slicedSurfaceFields.H"
#include "wedgeFvPatch.H"
#include "syncTools.H"

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
            mesh.Vsc0()().field()*rho.oldTime().field()
           *psi0*rDeltaT/mesh.Vsc()().field()
          + Su.field()
          - psiIf
        )/(rho.field()*rDeltaT - Sp.field());
    }
    else
    {
        psiIf =
        (
            rho.oldTime().field()*psi0*rDeltaT
          + Su.field()
          - psiIf
        )/(rho.field()*rDeltaT - Sp.field());
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
    class SuType,
    class PsiMaxType,
    class PsiMinType
>
void Foam::MULES::limiter
(
    scalarField& allLambda,
    const RdeltaTType& rDeltaT,
    const RhoType& rho,
    const volScalarField& psi,
    const surfaceScalarField& phiBD,
    const surfaceScalarField& phiCorr,
    const SpType& Sp,
    const SuType& Su,
    const PsiMaxType& psiMax,
    const PsiMinType& psiMin
)
{
    const scalarField& psiIf = psi;
    const volScalarField::Boundary& psiBf = psi.boundaryField();

    const fvMesh& mesh = psi.mesh();

    const dictionary& MULEScontrols = mesh.solverDict(psi.name());

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

    const scalarField& psi0 = psi.oldTime();

    const labelUList& owner = mesh.owner();
    const labelUList& neighb = mesh.neighbour();
    tmp<volScalarField::Internal> tVsc = mesh.Vsc();
    const scalarField& V = tVsc();

    const scalarField& phiBDIf = phiBD;
    const surfaceScalarField::Boundary& phiBDBf =
        phiBD.boundaryField();

    const scalarField& phiCorrIf = phiCorr;
    const surfaceScalarField::Boundary& phiCorrBf =
        phiCorr.boundaryField();

    slicedSurfaceScalarField lambda
    (
        IOobject
        (
            "lambda",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        mesh,
        dimless,
        allLambda,
        false   // Use slices for the couples
    );

    scalarField& lambdaIf = lambda;
    surfaceScalarField::Boundary& lambdaBf = lambda.boundaryFieldRef();

    scalarField psiMaxn(psiIf.size());
    scalarField psiMinn(psiIf.size());

    psiMaxn = psiMin;
    psiMinn = psiMax;

    scalarField sumPhiBD(psiIf.size(), 0.0);

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

        sumPhiBD[own] += phiBDIf[facei];
        sumPhiBD[nei] -= phiBDIf[facei];

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
        const scalarField& phiBDPf = phiBDBf[patchi];
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

            sumPhiBD[pfCelli] += phiBDPf[pFacei];

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

    psiMaxn = min(psiMaxn + extremaCoeff*(psiMax - psiMin), psiMax);
    psiMinn = max(psiMinn - extremaCoeff*(psiMax - psiMin), psiMin);

    if (smoothLimiter > small)
    {
        psiMaxn =
            min(smoothLimiter*psiIf + (1.0 - smoothLimiter)*psiMaxn, psiMax);
        psiMinn =
            max(smoothLimiter*psiIf + (1.0 - smoothLimiter)*psiMinn, psiMin);
    }

    if (mesh.moving())
    {
        tmp<volScalarField::Internal> V0 = mesh.Vsc0();

        psiMaxn =
            V
           *(
               (rho.field()*rDeltaT - Sp.field())*psiMaxn
             - Su.field()
            )
          - (V0().field()*rDeltaT)*rho.oldTime().field()*psi0
          + sumPhiBD;

        psiMinn =
            V
           *(
               Su.field()
             - (rho.field()*rDeltaT - Sp.field())*psiMinn
            )
          + (V0().field()*rDeltaT)*rho.oldTime().field()*psi0
          - sumPhiBD;
    }
    else
    {
        psiMaxn =
            V
           *(
               (rho.field()*rDeltaT - Sp.field())*psiMaxn
             - Su.field()
             - (rho.oldTime().field()*rDeltaT)*psi0
            )
          + sumPhiBD;

        psiMinn =
            V
           *(
               Su.field()
             - (rho.field()*rDeltaT - Sp.field())*psiMinn
             + (rho.oldTime().field()*rDeltaT)*psi0
            )
          - sumPhiBD;
    }

    scalarField sumlPhip(psiIf.size());
    scalarField mSumlPhim(psiIf.size());

    for (int j=0; j<nLimiterIter; j++)
    {
        sumlPhip = 0;
        mSumlPhim = 0;

        forAll(lambdaIf, facei)
        {
            const label own = owner[facei];
            const label nei = neighb[facei];

            scalar lambdaPhiCorrf = lambdaIf[facei]*phiCorrIf[facei];

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

            const labelList& pFaceCells = mesh.boundary()[patchi].faceCells();

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

        forAll(sumlPhip, celli)
        {
            sumlPhip[celli] =
                max(min
                (
                    (sumlPhip[celli] + psiMaxn[celli])
                   /(mSumPhim[celli] + rootVSmall),
                    1.0), 0.0
                );

            mSumlPhim[celli] =
                max(min
                (
                    (mSumlPhim[celli] + psiMinn[celli])
                   /(sumPhip[celli] + rootVSmall),
                    1.0), 0.0
                );
        }

        const scalarField& lambdam = sumlPhip;
        const scalarField& lambdap = mSumlPhim;

        forAll(lambdaIf, facei)
        {
            if (phiCorrIf[facei] > 0)
            {
                lambdaIf[facei] = min
                (
                    lambdaIf[facei],
                    min(lambdap[owner[facei]], lambdam[neighb[facei]])
                );
            }
            else
            {
                lambdaIf[facei] = min
                (
                    lambdaIf[facei],
                    min(lambdam[owner[facei]], lambdap[neighb[facei]])
                );
            }
        }

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

                forAll(lambdaPf, pFacei)
                {
                    const label pfCelli = pFaceCells[pFacei];

                    if (phiCorrfPf[pFacei] > 0)
                    {
                        lambdaPf[pFacei] =
                            min(lambdaPf[pFacei], lambdap[pfCelli]);
                    }
                    else
                    {
                        lambdaPf[pFacei] =
                            min(lambdaPf[pFacei], lambdam[pfCelli]);
                    }
                }
            }
        }

        syncTools::syncFaceList(mesh, allLambda, minEqOp<scalar>());
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

    surfaceScalarField::Boundary& phiBDBf = phiBD.boundaryFieldRef();
    const surfaceScalarField::Boundary& phiPsiBf = phiPsi.boundaryField();

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

    scalarField allLambda(mesh.nFaces(), 1.0);

    slicedSurfaceScalarField lambda
    (
        IOobject
        (
            "lambda",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        mesh,
        dimless,
        allLambda,
        false   // Use slices for the couples
    );

    limiter
    (
        allLambda,
        rDeltaT,
        rho,
        psi,
        phiBD,
        phiCorr,
        Sp,
        Su,
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
