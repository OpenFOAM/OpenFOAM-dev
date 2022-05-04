/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2022 OpenFOAM Foundation
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

#include "fvcMeshPhi.H"
#include "fvMesh.H"
#include "ddtScheme.H"
#include "surfaceInterpolate.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::tmp<Foam::surfaceScalarField> Foam::fvc::meshPhi
(
    const volVectorField& vf
)
{
    return fv::ddtScheme<vector>::New
    (
        vf.mesh(),
        vf.mesh().schemes().ddt("ddt(" + vf.name() + ')')
    ).ref().meshPhi(vf);
}


Foam::tmp<Foam::scalarField> Foam::fvc::meshPhi
(
    const volVectorField& vf,
    const label patchi
)
{
    return fv::ddtScheme<vector>::New
    (
        vf.mesh(),
        vf.mesh().schemes().ddt("ddt(" + vf.name() + ')')
    ).ref().meshPhi(vf, patchi);
}


Foam::tmp<Foam::surfaceScalarField> Foam::fvc::meshPhi
(
    const dimensionedScalar& rho,
    const volVectorField& vf
)
{
    return fv::ddtScheme<vector>::New
    (
        vf.mesh(),
        vf.mesh().schemes().ddt("ddt(" + rho.name() + ',' + vf.name() + ')')
    ).ref().meshPhi(vf);
}


Foam::tmp<Foam::surfaceScalarField> Foam::fvc::meshPhi
(
    const volScalarField& rho,
    const volVectorField& vf
)
{
    return fv::ddtScheme<vector>::New
    (
        vf.mesh(),
        vf.mesh().schemes().ddt("ddt(" + rho.name() + ',' + vf.name() + ')')
    ).ref().meshPhi(vf);
}


void Foam::fvc::makeRelative
(
    surfaceScalarField& phi,
    const volVectorField& U
)
{
    if (phi.mesh().moving())
    {
        phi -= fvc::meshPhi(U);
    }
}

void Foam::fvc::makeRelative
(
    surfaceScalarField& phi,
    const dimensionedScalar& rho,
    const volVectorField& U
)
{
    if (phi.mesh().moving())
    {
        phi -= rho*fvc::meshPhi(rho, U);
    }
}

void Foam::fvc::makeRelative
(
    surfaceScalarField& phi,
    const volScalarField& rho,
    const volVectorField& U
)
{
    if (phi.mesh().moving())
    {
        phi -= fvc::interpolate(rho)*fvc::meshPhi(rho, U);
    }
}


void Foam::fvc::makeAbsolute
(
    surfaceScalarField& phi,
    const volVectorField& U
)
{
    if (phi.mesh().moving())
    {
        phi += fvc::meshPhi(U);
    }
}

void Foam::fvc::makeAbsolute
(
    surfaceScalarField& phi,
    const dimensionedScalar& rho,
    const volVectorField& U
)
{
    if (phi.mesh().moving())
    {
        phi += rho*fvc::meshPhi(rho, U);
    }
}

void Foam::fvc::makeAbsolute
(
    surfaceScalarField& phi,
    const volScalarField& rho,
    const volVectorField& U
)
{
    if (phi.mesh().moving())
    {
        phi += fvc::interpolate(rho)*fvc::meshPhi(rho, U);
    }
}


Foam::tmp<Foam::surfaceScalarField> Foam::fvc::relative
(
    const tmp<surfaceScalarField>& tphi,
    const volVectorField& U
)
{
    if (tphi().mesh().moving())
    {
        return tphi - fvc::meshPhi(U);
    }
    else
    {
        return tmp<surfaceScalarField>(tphi, true);
    }
}


Foam::tmp<Foam::surfaceScalarField> Foam::fvc::relative
(
    const tmp<surfaceScalarField>& tphi,
    const volScalarField& rho,
    const volVectorField& U
)
{
    if (tphi().mesh().moving())
    {
        return tphi - fvc::interpolate(rho)*fvc::meshPhi(rho, U);
    }
    else
    {
        return tmp<surfaceScalarField>(tphi, true);
    }
}


Foam::tmp<Foam::surfaceScalarField> Foam::fvc::absolute
(
    const tmp<surfaceScalarField>& tphi,
    const volVectorField& U
)
{
    if (tphi().mesh().moving())
    {
        const word phiName(tphi().name());

        return surfaceScalarField::New
        (
            phiName,
            tphi + fvc::meshPhi(U)
        );
    }
    else
    {
        return tmp<surfaceScalarField>(tphi, true);
    }
}


Foam::tmp<Foam::surfaceScalarField> Foam::fvc::absolute
(
    const tmp<surfaceScalarField>& tphi,
    const volScalarField& rho,
    const volVectorField& U
)
{
    if (tphi().mesh().moving())
    {
        const word phiName(tphi().name());

        return surfaceScalarField::New
        (
            phiName,
            tphi + fvc::interpolate(rho)*fvc::meshPhi(rho, U)
        );
    }
    else
    {
        return tmp<surfaceScalarField>(tphi, true);
    }
}


Foam::tmp<Foam::surfaceScalarField> Foam::fvc::absolute
(
    const tmp<surfaceScalarField>& tphi,
    const volScalarField& alpha,
    const volScalarField& rho,
    const volVectorField& U
)
{
    if (tphi().mesh().moving())
    {
        const word phiName(tphi().name());

        return surfaceScalarField::New
        (
            phiName,
            tphi
          + fvc::interpolate(alpha)*fvc::interpolate(rho)*fvc::meshPhi(rho, U)
        );
    }
    else
    {
        return tmp<surfaceScalarField>(tphi, true);
    }
}


void Foam::fvc::correctUf
(
    autoPtr<surfaceVectorField>& Uf,
    const volVectorField& U,
    const surfaceScalarField& phi
)
{
    correctUf(Uf, U, phi, NullMRF());
}


void Foam::fvc::correctRhoUf
(
    autoPtr<surfaceVectorField>& rhoUf,
    const volScalarField& rho,
    const volVectorField& U,
    const surfaceScalarField& phi
)
{
    correctRhoUf(rhoUf, rho, U, phi, NullMRF());
}


// ************************************************************************* //
