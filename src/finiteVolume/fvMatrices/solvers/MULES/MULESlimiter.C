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
#include "volFields.H"
#include "surfaceFields.H"
#include "wedgeFvPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

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
    const control& controls,
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

    const scalar boundaryDeltaExtremaCoeff
    (
        max(controls.boundaryExtremaCoeff - controls.extremaCoeff, 0)
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

    scalarField phiCorrNorm;
    if (controls.tol != 0)
    {
        phiCorrNorm = (V*(rho.primitiveField()*rDeltaT - Sp.primitiveField()));
    }

    scalarField sumPhip(psiIf.size(), 0.0);
    scalarField mSumPhim(psiIf.size(), 0.0);

    if (controls.globalBounds)
    {
        psiMaxn = psiMax;
        psiMinn = psiMin;

        forAll(phiCorrIf, facei)
        {
            const label own = owner[facei];
            const label nei = neighb[facei];

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
            const scalarField& phiCorrPf = phiCorrBf[patchi];
            const labelList& pFaceCells = mesh.boundary()[patchi].faceCells();

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
    }
    else
    {
        psiMaxn = psiMin;
        psiMinn = psiMax;

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

        if (controls.extremaCoeff > 0)
        {
            psiMaxn = min
            (
                psiMaxn + controls.extremaCoeff*(psiMax - psiMin),
                psiMax
            );

            psiMinn = max
            (
                psiMinn - controls.extremaCoeff*(psiMax - psiMin),
                psiMin
            );
        }
        else
        {
            psiMaxn = min(psiMaxn, psiMax);
            psiMinn = max(psiMinn, psiMin);
        }

        if (controls.smoothingCoeff > small)
        {
            psiMaxn = min
            (
                controls.smoothingCoeff*psiIf
              + (1.0 - controls.smoothingCoeff)*psiMaxn,
                psiMax
            );

            psiMinn = max
            (
                controls.smoothingCoeff*psiIf
              + (1.0 - controls.smoothingCoeff)*psiMinn,
                psiMin
            );
        }
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
    if (controls.tol != 0)
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

    for (int j=0; j<controls.nIter; j++)
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

            if (controls.tol > 0)
            {
                const scalar phiCorrRes =
                    mag(phiCorrIf[facei])
                   /min(phiCorrNorm[owner[facei]], phiCorrNorm[neighb[facei]]);

                if (phiCorrRes > controls.tol)
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

                if (controls.tol > 0)
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

                        if (controls.tol > 0)
                        {
                            const scalar phiCorrRes =
                                mag(phiCorrfPf[pFacei])/phiCorrNorm[pfCelli];

                            if (phiCorrRes > controls.tol)
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

                if (controls.tol > 0)
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

                        if (phiCorrRes > controls.tol)
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
        if (controls.tol != 0)
        {
            reduce(maxDeltaLambdaPhiCorrRes, maxOp<scalar>());

            if (debug)
            {
                Info<< "MULES: maxDeltaLambdaPhiCorrRes "
                    << maxDeltaLambdaPhiCorrRes << endl;
            }

            if (maxDeltaLambdaPhiCorrRes < controls.tol) break;
        }
    }
}


// ************************************************************************* //
