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

template<class RdeltaTType, class RhoType, class SpType>
void Foam::MULES::correct
(
    const RdeltaTType& rDeltaT,
    const RhoType& rho,
    volScalarField& psi,
    const surfaceScalarField& phiCorr,
    const SpType& Sp
)
{
    Info<< "MULES: Correcting " << psi.name() << endl;

    scalarField psiIf(psi.size(), 0);
    fvc::surfaceIntegrate(psiIf, phiCorr);

    psi.primitiveFieldRef() =
    (
        (rho.primitiveField()*rDeltaT - Sp.primitiveField())
       *psi.primitiveField()
      - psiIf
    )/(rho.primitiveField()*rDeltaT - Sp.primitiveField());

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
    correct(rho, psi, phiCorr, zeroField());
}


template<class RhoType, class SpType>
void Foam::MULES::correct
(
    const RhoType& rho,
    volScalarField& psi,
    const surfaceScalarField& phiCorr,
    const SpType& Sp
)
{
    const fvMesh& mesh = psi.mesh();

    if (fv::localEulerDdt::enabled(mesh))
    {
        const volScalarField& rDeltaT = fv::localEulerDdt::localRDeltaT(mesh);
        correct(rDeltaT, rho, psi, phiCorr, Sp);
    }
    else
    {
        const scalar rDeltaT = 1.0/mesh.time().deltaTValue();
        correct(rDeltaT, rho, psi, phiCorr, Sp);
    }
}


template<class RhoType, class PsiMaxType, class PsiMinType>
void Foam::MULES::correct
(
    const control& controls,
    const RhoType& rho,
    volScalarField& psi,
    const surfaceScalarField& phiBD,
    surfaceScalarField& phiCorr,
    const PsiMaxType& psiMax,
    const PsiMinType& psiMin
)
{
    correct
    (
        controls,
        rho,
        psi,
        phiBD,
        phiCorr,
        zeroField(),
        psiMax,
        psiMin
    );
}


template
<
    class RhoType,
    class SpType,
    class PsiMaxType,
    class PsiMinType
>
void Foam::MULES::correct
(
    const control& controls,
    const RhoType& rho,
    volScalarField& psi,
    const surfaceScalarField& phiBD,
    surfaceScalarField& phiCorr,
    const SpType& Sp,
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
            controls,
            rDeltaT,
            rho,
            psi,
            phiBD,
            phiCorr,
            Sp,
            psiMax,
            psiMin
        );

        correct(rDeltaT, rho, psi, phiCorr, Sp);
    }
    else
    {
        const scalar rDeltaT = 1.0/mesh.time().deltaTValue();

        limitCorr
        (
            controls,
            rDeltaT,
            rho,
            psi,
            phiBD,
            phiCorr,
            Sp,
            psiMax,
            psiMin
        );

        correct(rDeltaT, rho, psi, phiCorr, Sp);
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
void Foam::MULES::limitCorr
(
    const control& controls,
    const RdeltaTType& rDeltaT,
    const RhoType& rho,
    const volScalarField& psi,
    const surfaceScalarField& phiBD,
    surfaceScalarField& phiCorr,
    const SpType& Sp,
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
       *(rho.primitiveField()*rDeltaT - Sp.primitiveField())
       *psi.primitiveField()

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

    phiCorr *= lambda;
}


template
<
    class RhoType,
    class SpType,
    class PsiMaxType,
    class PsiMinType
>
void Foam::MULES::limitCorr
(
    const control& controls,
    const RhoType& rho,
    const volScalarField& psi,
    const surfaceScalarField& phiBD,
    surfaceScalarField& phiCorr,
    const SpType& Sp,
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
            controls,
            rDeltaT,
            rho,
            psi,
            phiBD,
            phiCorr,
            Sp,
            psiMax,
            psiMin
        );
    }
    else
    {
        const scalar rDeltaT = 1.0/mesh.time().deltaTValue();

        limitCorr
        (
            controls,
            rDeltaT,
            rho,
            psi,
            phiBD,
            phiCorr,
            Sp,
            psiMax,
            psiMin
        );
    }
}


// ************************************************************************* //
