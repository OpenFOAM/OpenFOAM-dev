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
#include "linear.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class RdeltaTType, class RhoType, class SpType, class SuType>
void Foam::MULES::explicitSolve
(
    const RdeltaTType& rDeltaT,
    const RhoType& rho,
    volScalarField& psi,
    const surfaceScalarField& psiPhi,
    const SpType& Sp,
    const SuType& Su
)
{
    Info<< "MULES: Solving for " << psi.name() << endl;

    const fvMesh& mesh = psi.mesh();

    scalarField& psiIf = psi;
    const scalarField& psi0 = psi.oldTime();

    psiIf = 0.0;
    fvc::surfaceIntegrate(psiIf, psiPhi);

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
    const surfaceScalarField& psiPhi
)
{
    explicitSolve(rho, psi, psiPhi, zeroField(), zeroField());
}


template<class RhoType, class SpType, class SuType>
void Foam::MULES::explicitSolve
(
    const RhoType& rho,
    volScalarField& psi,
    const surfaceScalarField& psiPhi,
    const SpType& Sp,
    const SuType& Su
)
{
    const fvMesh& mesh = psi.mesh();

    if (fv::localEulerDdt::enabled(mesh))
    {
        const volScalarField& rDeltaT = fv::localEulerDdt::localRDeltaT(mesh);
        explicitSolve(rDeltaT, rho, psi, psiPhi, Sp, Su);
    }
    else
    {
        const scalar rDeltaT = 1.0/mesh.time().deltaTValue();
        explicitSolve(rDeltaT, rho, psi, psiPhi, Sp, Su);
    }
}


template<class RhoType, class PsiMaxType, class PsiMinType>
void Foam::MULES::explicitSolve
(
    const control& controls,
    const RhoType& rho,
    volScalarField& psi,
    const surfaceScalarField& phiBD,
    surfaceScalarField& psiPhi,
    const PsiMaxType& psiMax,
    const PsiMinType& psiMin
)
{
    explicitSolve
    (
        controls,
        rho,
        psi,
        phiBD,
        psiPhi,
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
    const control& controls,
    const RhoType& rho,
    volScalarField& psi,
    const surfaceScalarField& phi,
    surfaceScalarField& psiPhi,
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
        limit
        (
            controls,
            rDeltaT,
            rho,
            psi,
            phi,
            psiPhi,
            Sp,
            Su,
            psiMax,
            psiMin,
            false
        );
        explicitSolve(rDeltaT, rho, psi, psiPhi, Sp, Su);
    }
    else
    {
        const scalar rDeltaT = 1.0/mesh.time().deltaTValue();
        limit
        (
            controls,
            rDeltaT,
            rho,
            psi,
            phi,
            psiPhi,
            Sp,
            Su,
            psiMax,
            psiMin,
            false
        );
        explicitSolve(rDeltaT, rho, psi, psiPhi, Sp, Su);
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
    const control& controls,
    const RdeltaTType& rDeltaT,
    const RhoType& rho,
    const volScalarField& psi,
    const surfaceScalarField& phi,
    surfaceScalarField& psiPhi,
    const SpType& Sp,
    const SuType& Su,
    const PsiMaxType& psiMax,
    const PsiMinType& psiMin,
    const bool returnCorr
)
{
    surfaceScalarField phiBD(upwind<scalar>(psi.mesh(), phi).flux(psi));

    surfaceScalarField::Boundary& phiBDBf = phiBD.boundaryFieldRef();
    const surfaceScalarField::Boundary& psiPhiBf = psiPhi.boundaryField();

    forAll(phiBDBf, patchi)
    {
        fvsPatchScalarField& phiBDPf = phiBDBf[patchi];

        if (!phiBDPf.coupled())
        {
            phiBDPf = psiPhiBf[patchi];
        }
    }

    limit
    (
        controls,
        rDeltaT,
        rho,
        psi,
        phi,
        phiBD,
        psiPhi,
        Sp,Su,
        psiMax,
        psiMin,
        returnCorr
    );
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
    const control& controls,
    const RdeltaTType& rDeltaT,
    const RhoType& rho,
    const volScalarField& psi,
    const surfaceScalarField& phi,
    const surfaceScalarField& phiBD,
    surfaceScalarField& psiPhi,
    const SpType& Sp,
    const SuType& Su,
    const PsiMaxType& psiMax,
    const PsiMinType& psiMin,
    const bool returnCorr
)
{
    const fvMesh& mesh = psi.mesh();

    const labelUList& owner = mesh.owner();
    const labelUList& neighb = mesh.neighbour();

    const scalarField& phiBDIf = phiBD;
    const surfaceScalarField::Boundary& phiBDBf = phiBD.boundaryField();

    surfaceScalarField& phiCorr = psiPhi;
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
        controls,
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
        psiPhi = phiBD + lambda*phiCorr;
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
    const control& controls,
    const RhoType& rho,
    const volScalarField& psi,
    const surfaceScalarField& phi,
    surfaceScalarField& psiPhi,
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
        limit
        (
            controls,
            rDeltaT,
            rho,
            psi,
            phi,
            psiPhi,
            Sp,
            Su,
            psiMax,
            psiMin,
            rtnCorr
        );
    }
    else
    {
        const scalar rDeltaT = 1.0/mesh.time().deltaTValue();
        limit
        (
            controls,
            rDeltaT,
            rho,
            psi,
            phi,
            psiPhi,
            Sp,
            Su,
            psiMax,
            psiMin,
            rtnCorr
        );
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
    const control& controls,
    const RhoType& rho,
    const volScalarField& psi,
    const surfaceScalarField& phi,
    const surfaceScalarField& phiBD,
    surfaceScalarField& psiPhi,
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
        limit
        (
            controls,
            rDeltaT,
            rho,
            psi,
            phi,
            phiBD,
            psiPhi,
            Sp,
            Su,
            psiMax,
            psiMin,
            rtnCorr
        );
    }
    else
    {
        const scalar rDeltaT = 1.0/mesh.time().deltaTValue();
        limit
        (
            controls,
            rDeltaT,
            rho,
            psi,
            phi,
            phiBD,
            psiPhi,
            Sp,
            Su,
            psiMax,
            psiMin,
            rtnCorr
        );
    }
}


template
<
    template<class> class ConstraintList,
    class ConstraintField,
    template<class> class PsiPhiList
>
void Foam::MULES::limitSum
(
    const ConstraintList<const ConstraintField>& constraints,
    PsiPhiList<surfaceScalarField>& psiPhiCorrs,
    const surfaceScalarField& phi
)
{
    {
        UPtrList<scalarField> psiPhiCorrsInternal(psiPhiCorrs.size());

        forAll(psiPhiCorrsInternal, phasei)
        {
            psiPhiCorrsInternal.set(phasei, &psiPhiCorrs[phasei]);
        }

        limitSum(psiPhiCorrsInternal);
    }

    const surfaceScalarField::Boundary& phibf = phi.boundaryField();

    forAll(phibf, patchi)
    {
        label nFixed = 0;

        forAll(constraints, phasei)
        {
            if (constraints[phasei].boundaryField()[patchi].fixesValue())
            {
                nFixed++;
            }
        }

        if (nFixed == 0)
        {
            UPtrList<scalarField> psiPhiCorrsPatch(psiPhiCorrs.size());

            forAll(psiPhiCorrsPatch, phasei)
            {
                psiPhiCorrsPatch.set
                (
                    phasei,
                    &psiPhiCorrs[phasei].boundaryFieldRef()[patchi]
                );
            }

            limitSum(psiPhiCorrsPatch);
        }
        else if (nFixed < psiPhiCorrs.size())
        {
            UPtrList<scalarField> psiPhiCorrsPatch(psiPhiCorrs.size());
            UPtrList<scalarField> psiPhiCorrsFixedPatch(psiPhiCorrs.size());

            label i = 0;
            label fixedi = 0;

            forAll(psiPhiCorrsPatch, phasei)
            {
                if (constraints[phasei].boundaryField()[patchi].fixesValue())
                {
                    psiPhiCorrsFixedPatch.set
                    (
                        fixedi++,
                        &psiPhiCorrs[phasei].boundaryFieldRef()[patchi]
                    );
                }
                else
                {
                    psiPhiCorrsPatch.set
                    (
                        i++,
                        &psiPhiCorrs[phasei].boundaryFieldRef()[patchi]
                    );
                }
            }

            psiPhiCorrsPatch.setSize(i);
            psiPhiCorrsFixedPatch.setSize(fixedi);

            limitSum(psiPhiCorrsPatch, psiPhiCorrsFixedPatch);
        }
    }
}


template
<
    template<class> class ConstraintList,
    class ConstraintField,
    template<class> class PsiPhiList
>
void Foam::MULES::limitSum
(
    const PtrList<surfaceScalarField>& alphaPhiBDs,
    const ConstraintList<const ConstraintField>& constraints,
    PsiPhiList<surfaceScalarField>& psiPhis,
    const surfaceScalarField& phi
)
{
    forAll(psiPhis, phasei)
    {
        psiPhis[phasei] -= alphaPhiBDs[phasei];
    }

    limitSum(constraints, psiPhis, phi);

    forAll(psiPhis, phasei)
    {
        psiPhis[phasei] += alphaPhiBDs[phasei];
    }
}


template
<
    template<class> class AlphaList,
    template<class> class ConstraintList,
    class ConstraintField,
    template<class> class PsiPhiList
>
void Foam::MULES::limitSum
(
    const AlphaList<const volScalarField>& alphas,
    const ConstraintList<const ConstraintField>& constraints,
    PsiPhiList<surfaceScalarField>& psiPhis,
    const surfaceScalarField& phi
)
{
    PtrList<surfaceScalarField> alphaPhiBDs(psiPhis.size());

    forAll(psiPhis, phasei)
    {
        alphaPhiBDs.set
        (
            phasei,
            upwind<scalar>(phi.mesh(), phi).flux(alphas[phasei])
        );
    }

    limitSum(alphaPhiBDs, constraints, psiPhis, phi);
}


// ************************************************************************* //
