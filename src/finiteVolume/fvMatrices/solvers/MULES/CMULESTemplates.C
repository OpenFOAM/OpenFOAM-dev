/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2013-2025 OpenFOAM Foundation
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
            rho.primitiveField()*psi.primitiveField()*rDeltaT
          + Su.primitiveField()
          - psiIf
        )/(rho.primitiveField()*rDeltaT - Sp.primitiveField());
    }
    else
    {
        psi.primitiveFieldRef() =
        (
            rho.primitiveField()*psi.primitiveField()*rDeltaT
          + Su.primitiveField()
          - psiIf
        )/(rho.primitiveField()*rDeltaT - Sp.primitiveField());
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
    const surfaceScalarField& phiBD,
    surfaceScalarField& phiCorr,
    const PsiMaxType& psiMax,
    const PsiMinType& psiMin
)
{
    correct(rho, psi, phiBD, phiCorr, zeroField(), zeroField(), psiMax, psiMin);
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
    const surfaceScalarField& phiBD,
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
            phiBD,
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
            phiBD,
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
void Foam::MULES::limitCorr
(
    const RdeltaTType& rDeltaT,
    const RhoType& rho,
    const volScalarField& psi,
    const surfaceScalarField& phiBD,
    surfaceScalarField& phiCorr,
    const SpType& Sp,
    const SuType& Su,
    const PsiMaxType& psiMax,
    const PsiMinType& psiMin
)
{
    const fvMesh& mesh = psi.mesh();

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

    // Correction equation source
    scalarField SuCorr
    (
        mesh.Vsc()().primitiveField()
       *(
            (rho.primitiveField()*rDeltaT)*psi.primitiveField()
          + Su.primitiveField()
        )
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

    phiCorr *= lambda;
}


// ************************************************************************* //
