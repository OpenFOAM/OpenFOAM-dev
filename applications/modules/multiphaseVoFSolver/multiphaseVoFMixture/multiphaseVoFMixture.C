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

#include "multiphaseVoFMixture.H"
#include "correctContactAngle.H"
#include "surfaceInterpolate.H"
#include "fvcGrad.H"
#include "fvcSnGrad.H"
#include "fvcDiv.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(multiphaseVoFMixture, 0);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::multiphaseVoFMixture::calcAlphas()
{
    scalar level = 0.0;
    alphas_ == 0.0;

    forAll(phases_, phasei)
    {
        alphas_ += level*phases_[phasei];
        level += 1.0;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::multiphaseVoFMixture::multiphaseVoFMixture
(
    const fvMesh& mesh,
    const VoFphase::iNew& inew
)
:
    VoFMixture(mesh),
    IOdictionary
    (
        IOobject
        (
            "phaseProperties",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    ),

    phases_(lookup("phases"), inew),

    mesh_(mesh),

    alphas_
    (
        IOobject
        (
            "alphas",
            mesh_.time().name(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimless, 0)
    ),

    sigmas_(lookup("sigmas")),
    dimSigma_(1, 0, -2, 0, 0),

    deltaN_
    (
        "deltaN",
        1e-8/pow(average(mesh_.V()), 1.0/3.0)
    )
{
    calcAlphas();
    alphas_.write();
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::surfaceScalarField>
Foam::multiphaseVoFMixture::surfaceTensionForce(const volVectorField& U) const
{
    tmp<surfaceScalarField> tstf
    (
        surfaceScalarField::New
        (
            "surfaceTensionForce",
            mesh_,
            dimensionedScalar(dimensionSet(1, -2, -2, 0, 0), 0)
        )
    );

    surfaceScalarField& stf = tstf.ref();

    forAll(phases_, phasei)
    {
        const VoFphase& alpha1 = phases_[phasei];

        for (label phasej = phasei+1; phasej<phases_.size(); phasej++)
        {
            const VoFphase& alpha2 = phases_[phasej];

            sigmaTable::const_iterator sigma =
                sigmas_.find(interfacePair(alpha1, alpha2));

            if (sigma == sigmas_.end())
            {
                FatalErrorInFunction
                    << "Cannot find interface " << interfacePair(alpha1, alpha2)
                    << " in list of sigma values"
                    << exit(FatalError);
            }

            stf += dimensionedScalar(dimSigma_, sigma())
               *fvc::interpolate(K(alpha1, alpha2, U))*
                (
                    fvc::interpolate(alpha2)*fvc::snGrad(alpha1)
                  - fvc::interpolate(alpha1)*fvc::snGrad(alpha2)
                );
        }
    }

    return tstf;
}


Foam::tmp<Foam::surfaceVectorField> Foam::multiphaseVoFMixture::nHatfv
(
    const volScalarField& alpha1,
    const volScalarField& alpha2
) const
{
    /*
    // Cell gradient of alpha
    volVectorField gradAlpha =
        alpha2*fvc::grad(alpha1) - alpha1*fvc::grad(alpha2);

    // Interpolated face-gradient of alpha
    surfaceVectorField gradAlphaf = fvc::interpolate(gradAlpha);
    */

    surfaceVectorField gradAlphaf
    (
        fvc::interpolate(alpha2)*fvc::interpolate(fvc::grad(alpha1))
      - fvc::interpolate(alpha1)*fvc::interpolate(fvc::grad(alpha2))
    );

    // Face unit interface normal
    return gradAlphaf/(mag(gradAlphaf) + deltaN_);
}


Foam::tmp<Foam::surfaceScalarField> Foam::multiphaseVoFMixture::nHatf
(
    const volScalarField& alpha1,
    const volScalarField& alpha2
) const
{
    // Face unit interface normal flux
    return nHatfv(alpha1, alpha2) & mesh_.Sf();
}


Foam::tmp<Foam::volScalarField> Foam::multiphaseVoFMixture::K
(
    const VoFphase& alpha1,
    const VoFphase& alpha2,
    const volVectorField& U
) const
{
    tmp<surfaceVectorField> tnHatfv = nHatfv(alpha1, alpha2);

    correctContactAngle
    (
        alpha1,
        alpha2,
        U.boundaryField(),
        deltaN_,
        tnHatfv.ref().boundaryFieldRef()
    );

    // Simple expression for curvature
    return -fvc::div(tnHatfv & mesh_.Sf());
}


Foam::tmp<Foam::volScalarField>
Foam::multiphaseVoFMixture::nearInterface() const
{
    tmp<volScalarField> tnearInt
    (
        volScalarField::New
        (
            "nearInterface",
            mesh_,
            dimensionedScalar(dimless, 0)
        )
    );

    forAll(phases_, phasei)
    {
        tnearInt.ref() = max
        (
            tnearInt(),
            pos0(phases_[phasei] - 0.01)*pos0(0.99 - phases_[phasei])
        );
    }

    return tnearInt;
}


bool Foam::multiphaseVoFMixture::read()
{
    if (regIOobject::read())
    {
        lookup("sigmas") >> sigmas_;
        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
