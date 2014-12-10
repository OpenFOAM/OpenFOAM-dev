/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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
    const scalar rDeltaT = 1.0/mesh.time().deltaTValue();
    explicitSolve(rDeltaT, rho, psi, phiPsi, Sp, Su);
}


template<class RhoType, class SpType, class SuType>
void Foam::MULES::explicitSolve
(
    const RhoType& rho,
    volScalarField& psi,
    const surfaceScalarField& phi,
    surfaceScalarField& phiPsi,
    const SpType& Sp,
    const SuType& Su,
    const scalar psiMax,
    const scalar psiMin
)
{
    const fvMesh& mesh = psi.mesh();
    const scalar rDeltaT = 1.0/mesh.time().deltaTValue();
    psi.correctBoundaryConditions();
    limit(rDeltaT, rho, psi, phi, phiPsi, Sp, Su, psiMax, psiMin, 3, false);
    explicitSolve(rDeltaT, rho, psi, phiPsi, Sp, Su);
}


template<class RhoType, class SpType, class SuType>
void Foam::MULES::explicitLTSSolve
(
    const RhoType& rho,
    volScalarField& psi,
    const surfaceScalarField& phi,
    surfaceScalarField& phiPsi,
    const SpType& Sp,
    const SuType& Su,
    const scalar psiMax,
    const scalar psiMin
)
{
    const fvMesh& mesh = psi.mesh();

    const volScalarField& rDeltaT =
        mesh.objectRegistry::lookupObject<volScalarField>("rSubDeltaT");

    psi.correctBoundaryConditions();
    limit(rDeltaT, rho, psi, phi, phiPsi, Sp, Su, psiMax, psiMin, 3, false);
    explicitSolve(rDeltaT, rho, psi, phiPsi, Sp, Su);
}


template<class RdeltaTType, class RhoType, class SpType, class SuType>
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
    const scalar psiMax,
    const scalar psiMin,
    const label nLimiterIter
)
{
    const scalarField& psiIf = psi;
    const volScalarField::GeometricBoundaryField& psiBf = psi.boundaryField();

    const scalarField& psi0 = psi.oldTime();

    const fvMesh& mesh = psi.mesh();

    const labelUList& owner = mesh.owner();
    const labelUList& neighb = mesh.neighbour();
    tmp<volScalarField::DimensionedInternalField> tVsc = mesh.Vsc();
    const scalarField& V = tVsc();

    const scalarField& phiBDIf = phiBD;
    const surfaceScalarField::GeometricBoundaryField& phiBDBf =
        phiBD.boundaryField();

    const scalarField& phiCorrIf = phiCorr;
    const surfaceScalarField::GeometricBoundaryField& phiCorrBf =
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
    surfaceScalarField::GeometricBoundaryField& lambdaBf =
        lambda.boundaryField();

    scalarField psiMaxn(psiIf.size(), psiMin);
    scalarField psiMinn(psiIf.size(), psiMax);

    scalarField sumPhiBD(psiIf.size(), 0.0);

    scalarField sumPhip(psiIf.size(), VSMALL);
    scalarField mSumPhim(psiIf.size(), VSMALL);

    forAll(phiCorrIf, facei)
    {
        label own = owner[facei];
        label nei = neighb[facei];

        psiMaxn[own] = max(psiMaxn[own], psiIf[nei]);
        psiMinn[own] = min(psiMinn[own], psiIf[nei]);

        psiMaxn[nei] = max(psiMaxn[nei], psiIf[own]);
        psiMinn[nei] = min(psiMinn[nei], psiIf[own]);

        sumPhiBD[own] += phiBDIf[facei];
        sumPhiBD[nei] -= phiBDIf[facei];

        scalar phiCorrf = phiCorrIf[facei];

        if (phiCorrf > 0.0)
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
                label pfCelli = pFaceCells[pFacei];

                psiMaxn[pfCelli] = max(psiMaxn[pfCelli], psiPNf[pFacei]);
                psiMinn[pfCelli] = min(psiMinn[pfCelli], psiPNf[pFacei]);
            }
        }
        else
        {
            forAll(phiCorrPf, pFacei)
            {
                label pfCelli = pFaceCells[pFacei];

                psiMaxn[pfCelli] = max(psiMaxn[pfCelli], psiPf[pFacei]);
                psiMinn[pfCelli] = min(psiMinn[pfCelli], psiPf[pFacei]);
            }
        }

        forAll(phiCorrPf, pFacei)
        {
            label pfCelli = pFaceCells[pFacei];

            sumPhiBD[pfCelli] += phiBDPf[pFacei];

            scalar phiCorrf = phiCorrPf[pFacei];

            if (phiCorrf > 0.0)
            {
                sumPhip[pfCelli] += phiCorrf;
            }
            else
            {
                mSumPhim[pfCelli] -= phiCorrf;
            }
        }
    }

    psiMaxn = min(psiMaxn, psiMax);
    psiMinn = max(psiMinn, psiMin);

    //scalar smooth = 0.5;
    //psiMaxn = min((1.0 - smooth)*psiIf + smooth*psiMaxn, psiMax);
    //psiMinn = max((1.0 - smooth)*psiIf + smooth*psiMinn, psiMin);

    if (mesh.moving())
    {
        tmp<volScalarField::DimensionedInternalField> V0 = mesh.Vsc0();

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
        sumlPhip = 0.0;
        mSumlPhim = 0.0;

        forAll(lambdaIf, facei)
        {
            label own = owner[facei];
            label nei = neighb[facei];

            scalar lambdaPhiCorrf = lambdaIf[facei]*phiCorrIf[facei];

            if (lambdaPhiCorrf > 0.0)
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
                label pfCelli = pFaceCells[pFacei];

                scalar lambdaPhiCorrf = lambdaPf[pFacei]*phiCorrfPf[pFacei];

                if (lambdaPhiCorrf > 0.0)
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
                   /(mSumPhim[celli] - SMALL),
                    1.0), 0.0
                );

            mSumlPhim[celli] =
                max(min
                (
                    (mSumlPhim[celli] + psiMinn[celli])
                   /(sumPhip[celli] + SMALL),
                    1.0), 0.0
                );
        }

        const scalarField& lambdam = sumlPhip;
        const scalarField& lambdap = mSumlPhim;

        forAll(lambdaIf, facei)
        {
            if (phiCorrIf[facei] > 0.0)
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
                    label pfCelli = pFaceCells[pFacei];

                    if (phiCorrfPf[pFacei] > 0.0)
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
            else
            {
                const labelList& pFaceCells =
                    mesh.boundary()[patchi].faceCells();
                const scalarField& phiBDPf = phiBDBf[patchi];
                const scalarField& phiCorrPf = phiCorrBf[patchi];

                forAll(lambdaPf, pFacei)
                {
                    // Limit outlet faces only
                    if ((phiBDPf[pFacei] + phiCorrPf[pFacei]) > SMALL*SMALL)
                    {
                        label pfCelli = pFaceCells[pFacei];

                        if (phiCorrfPf[pFacei] > 0.0)
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
        }

        syncTools::syncFaceList(mesh, allLambda, minEqOp<scalar>());
    }
}


template<class RdeltaTType, class RhoType, class SpType, class SuType>
void Foam::MULES::limit
(
    const RdeltaTType& rDeltaT,
    const RhoType& rho,
    const volScalarField& psi,
    const surfaceScalarField& phi,
    surfaceScalarField& phiPsi,
    const SpType& Sp,
    const SuType& Su,
    const scalar psiMax,
    const scalar psiMin,
    const label nLimiterIter,
    const bool returnCorr
)
{
    const fvMesh& mesh = psi.mesh();

    surfaceScalarField phiBD(upwind<scalar>(psi.mesh(), phi).flux(psi));

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
        psiMin,
        nLimiterIter
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


template<class SurfaceScalarFieldList>
void Foam::MULES::limitSum(SurfaceScalarFieldList& phiPsiCorrs)
{
    {
        UPtrList<scalarField> phiPsiCorrsInternal(phiPsiCorrs.size());
        forAll(phiPsiCorrs, phasei)
        {
            phiPsiCorrsInternal.set(phasei, &phiPsiCorrs[phasei]);
        }

        limitSum(phiPsiCorrsInternal);
    }

    surfaceScalarField::GeometricBoundaryField& bfld =
        phiPsiCorrs[0].boundaryField();

    forAll(bfld, patchi)
    {
        if (bfld[patchi].coupled())
        {
            UPtrList<scalarField> phiPsiCorrsPatch(phiPsiCorrs.size());
            forAll(phiPsiCorrs, phasei)
            {
                phiPsiCorrsPatch.set
                (
                    phasei,
                    &phiPsiCorrs[phasei].boundaryField()[patchi]
                );
            }

            limitSum(phiPsiCorrsPatch);
        }
    }
}


// ************************************************************************* //
