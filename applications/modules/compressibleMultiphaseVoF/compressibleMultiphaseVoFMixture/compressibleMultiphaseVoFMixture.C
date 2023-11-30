/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2023 OpenFOAM Foundation
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

#include "compressibleMultiphaseVoFMixture.H"
#include "surfaceInterpolate.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(compressibleMultiphaseVoFMixture, 0);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::compressibleMultiphaseVoFMixture::compressibleMultiphaseVoFMixture
(
    const fvMesh& mesh
)
:
    compressibleMultiphaseVoFMixtureThermo(mesh),

    multiphaseVoFMixture(mesh, compressibleVoFphase::iNew(mesh, T())),

    phases_(multiphaseVoFMixture::phases().convert<compressibleVoFphase>()),

    rho_
    (
        IOobject
        (
            "rho",
            mesh.time().name(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("rho", dimDensity, 0)
    )
{
    correct();
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::compressibleMultiphaseVoFMixture::incompressible() const
{
    bool incompressible = true;

    forAll(phases_, phasei)
    {
        incompressible =
            incompressible && phases_[phasei].thermo().incompressible();
    }

    return incompressible;
}


Foam::tmp<Foam::volScalarField>
Foam::compressibleMultiphaseVoFMixture::nu() const
{
    volScalarField mu(phases_[0].Alpha()*phases_[0].thermo().mu());

    for (label phasei=1; phasei<phases_.size(); phasei++)
    {
        mu += phases_[phasei].Alpha()*phases_[phasei].thermo().mu();
    }

    return mu/rho_;
}


Foam::tmp<Foam::scalarField> Foam::compressibleMultiphaseVoFMixture::nu
(
    const label patchi
) const
{
    scalarField mu
    (
        phases_[0].Alpha().boundaryField()[patchi]
       *phases_[0].thermo().mu().boundaryField()[patchi]
    );

    for (label phasei=1; phasei<phases_.size(); phasei++)
    {
        mu +=
            phases_[phasei].Alpha().boundaryField()[patchi]
           *phases_[phasei].thermo().mu().boundaryField()[patchi];
    }

    return mu/rho_.boundaryField()[patchi];
}


Foam::tmp<Foam::volScalarField>
Foam::compressibleMultiphaseVoFMixture::psiByRho() const
{
    tmp<volScalarField> tpsiByRho
    (
        phases_[0]*phases_[0].thermo().psi()/phases_[0].thermo().rho()
    );

    for (label phasei=1; phasei<phases_.size(); phasei++)
    {
        tpsiByRho.ref() +=
            phases_[phasei]*phases_[phasei].thermo().psi()
           /phases_[phasei].thermo().rho();
    }

    return tpsiByRho;
}


Foam::tmp<Foam::volScalarField> Foam::compressibleMultiphaseVoFMixture::alphaEff
(
    const volScalarField& nut
) const
{
    tmp<volScalarField> talphaEff
    (
        phases_[0]
       *(
           phases_[0].thermo().kappa()
         + phases_[0].thermo().rho()*phases_[0].thermo().Cp()*nut
        )/phases_[0].thermo().Cv()
    );

    for (label phasei=1; phasei<phases_.size(); phasei++)
    {
        talphaEff.ref() +=
            phases_[phasei]
           *(
               phases_[phasei].thermo().kappa()
             + phases_[phasei].thermo().rho()*phases_[phasei].thermo().Cp()*nut
            )/phases_[phasei].thermo().Cv();
    }

    return talphaEff;
}


Foam::tmp<Foam::volScalarField>
Foam::compressibleMultiphaseVoFMixture::rCv() const
{
    tmp<volScalarField> trCv(phases_[0]/phases_[0].thermo().Cv());

    for (label phasei=1; phasei<phases_.size(); phasei++)
    {
        trCv.ref() += phases_[phasei]/phases_[phasei].thermo().Cv();
    }

    return trCv;
}


void Foam::compressibleMultiphaseVoFMixture::correctThermo()
{
    forAll(phases_, phasei)
    {
        phases_[phasei].correct(p(), T());
    }
}


void Foam::compressibleMultiphaseVoFMixture::correct()
{
    rho_ = phases_[0]*phases_[0].thermo().rho();

    for (label phasei=1; phasei<phases_.size(); phasei++)
    {
        rho_ += phases_[phasei]*phases_[phasei].thermo().rho();
    }

    forAll(phases_, phasei)
    {
        phases_[phasei].Alpha() =
            phases_[phasei]*phases_[phasei].thermo().rho()/rho_;
    }

    calcAlphas();
}


void Foam::compressibleMultiphaseVoFMixture::correctRho
(
    const volScalarField& dp
)
{
    forAll(phases_, phasei)
    {
        phases_[phasei].thermo().rho() += phases_[phasei].thermo().psi()*dp;
    }
}


// ************************************************************************* //
