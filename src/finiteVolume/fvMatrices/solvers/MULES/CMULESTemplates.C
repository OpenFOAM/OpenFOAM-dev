/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2013-2019 OpenFOAM Foundation
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

#include "CMULES.H"
#include "fvcSurfaceIntegrate.H"
#include "localEulerDdtScheme.H"
#include "slicedSurfaceFields.H"
#include "wedgeFvPatch.H"
#include "syncTools.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class RdeltaTType, class RhoType, class SpType, class SuType>
void Foam::MULES::correct
(
    const RdeltaTType& rDeltaT,
    const RhoType& rho,
    volScalarField& psi,
    const surfaceScalarField& phiCorr,
    const SpType& Sp,
    const SuType& Su
)
{
    Info<< "MULES: Correcting " << psi.name() << endl;

    const fvMesh& mesh = psi.mesh();

    scalarField psiIf(psi.size(), 0);
    fvc::surfaceIntegrate(psiIf, phiCorr);

    if (mesh.moving())
    {
        psi.primitiveFieldRef() =
        (
            rho.field()*psi.primitiveField()*rDeltaT
          + Su.field()
          - psiIf
        )/(rho.field()*rDeltaT - Sp.field());
    }
    else
    {
        psi.primitiveFieldRef() =
        (
            rho.field()*psi.primitiveField()*rDeltaT
          + Su.field()
          - psiIf
        )/(rho.field()*rDeltaT - Sp.field());
    }

    psi.correctBoundaryConditions();
}


template<class RhoType>
void Foam::MULES::correct
(
    const RhoType& rho,
    volScalarField& psi,
    const surfaceScalarField& phiCorr
)
{
    correct(rho, psi, phiCorr, zeroField(), zeroField());
}


template<class RhoType, class SpType, class SuType>
void Foam::MULES::correct
(
    const RhoType& rho,
    volScalarField& psi,
    const surfaceScalarField& phiCorr,
    const SpType& Sp,
    const SuType& Su
)
{
    const fvMesh& mesh = psi.mesh();

    if (fv::localEulerDdt::enabled(mesh))
    {
        const volScalarField& rDeltaT = fv::localEulerDdt::localRDeltaT(mesh);
        correct(rDeltaT, rho, psi, phiCorr, Sp, Su);
    }
    else
    {
        const scalar rDeltaT = 1.0/mesh.time().deltaTValue();
        correct(rDeltaT, rho, psi, phiCorr, Sp, Su);
    }
}


template<class RhoType, class PsiMaxType, class PsiMinType>
void Foam::MULES::correct
(
    const RhoType& rho,
    volScalarField& psi,
    const surfaceScalarField& phi,
    surfaceScalarField& phiCorr,
    const PsiMaxType& psiMax,
    const PsiMinType& psiMin
)
{
    correct(rho, psi, phi, phiCorr, zeroField(), zeroField(), psiMax, psiMin);
}


template
<
    class RhoType,
    class SpType,
    class SuType,
    class PsiMaxType,
    class PsiMinType
>
void Foam::MULES::correct
(
    const RhoType& rho,
    volScalarField& psi,
    const surfaceScalarField& phi,
    surfaceScalarField& phiCorr,
    const SpType& Sp,
    const SuType& Su,
    const PsiMaxType& psiMax,
    const PsiMinType& psiMin
)
{
    const fvMesh& mesh = psi.mesh();

    if (fv::localEulerDdt::enabled(mesh))
    {
        const volScalarField& rDeltaT = fv::localEulerDdt::localRDeltaT(mesh);

        limitCorr
        (
            rDeltaT,
            rho,
            psi,
            phi,
            phiCorr,
            Sp,
            Su,
            psiMax,
            psiMin
        );

        correct(rDeltaT, rho, psi, phiCorr, Sp, Su);
    }
    else
    {
        const scalar rDeltaT = 1.0/mesh.time().deltaTValue();

        limitCorr
        (
            rDeltaT,
            rho,
            psi,
            phi,
            phiCorr,
            Sp,
            Su,
            psiMax,
            psiMin
        );

        correct(rDeltaT, rho, psi, phiCorr, Sp, Su);
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
void Foam::MULES::limiterCorr
(
    scalarField& allLambda,
    const RdeltaTType& rDeltaT,
    const RhoType& rho,
    const volScalarField& psi,
    const surfaceScalarField& phi,
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
        MULEScontrols.lookup<label>("nLimiterIter")
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

    const labelUList& owner = mesh.owner();
    const labelUList& neighb = mesh.neighbour();
    tmp<volScalarField::Internal> tVsc = mesh.Vsc();
    const scalarField& V = tVsc();

    const surfaceScalarField::Boundary& phiBf =
        phi.boundaryField();

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
    surfaceScalarField::Boundary& lambdaBf =
        lambda.boundaryFieldRef();

    scalarField psiMaxn(psiIf.size());
    scalarField psiMinn(psiIf.size());

    psiMaxn = psiMin;
    psiMinn = psiMax;

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
                label pfCelli = pFaceCells[pFacei];

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

    psiMaxn = min(psiMaxn + extremaCoeff*(psiMax - psiMin), psiMax);
    psiMinn = max(psiMinn - extremaCoeff*(psiMax - psiMin), psiMin);

    if (smoothLimiter > small)
    {
        psiMaxn =
            min(smoothLimiter*psiIf + (1.0 - smoothLimiter)*psiMaxn, psiMax);
        psiMinn =
            max(smoothLimiter*psiIf + (1.0 - smoothLimiter)*psiMinn, psiMin);
    }

    psiMaxn =
        V
       *(
           (rho.field()*rDeltaT - Sp.field())*psiMaxn
         - Su.field()
         - rho.field()*psi.primitiveField()*rDeltaT
        );

    psiMinn =
        V
       *(
           Su.field()
         - (rho.field()*rDeltaT - Sp.field())*psiMinn
         + rho.field()*psi.primitiveField()*rDeltaT
        );

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

            const labelList& pFaceCells = mesh.boundary()[patchi].faceCells();

            forAll(lambdaPf, pFacei)
            {
                label pfCelli = pFaceCells[pFacei];

                scalar lambdaPhiCorrf = lambdaPf[pFacei]*phiCorrfPf[pFacei];

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
            else
            {
                const labelList& pFaceCells =
                    mesh.boundary()[patchi].faceCells();
                const scalarField& phiPf = phiBf[patchi];

                forAll(lambdaPf, pFacei)
                {
                    // Limit outlet faces only
                    if ((phiPf[pFacei] + phiCorrfPf[pFacei]) > small*small)
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
void Foam::MULES::limitCorr
(
    const RdeltaTType& rDeltaT,
    const RhoType& rho,
    const volScalarField& psi,
    const surfaceScalarField& phi,
    surfaceScalarField& phiCorr,
    const SpType& Sp,
    const SuType& Su,
    const PsiMaxType& psiMax,
    const PsiMinType& psiMin
)
{
    const fvMesh& mesh = psi.mesh();

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

    limiterCorr
    (
        allLambda,
        rDeltaT,
        rho,
        psi,
        phi,
        phiCorr,
        Sp,
        Su,
        psiMax,
        psiMin
    );

    phiCorr *= lambda;
}


// ************************************************************************* //
