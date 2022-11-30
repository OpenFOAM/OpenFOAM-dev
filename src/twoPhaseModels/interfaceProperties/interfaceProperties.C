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

#include "interfaceProperties.H"
#include "alphaContactAngleFvPatchScalarField.H"
#include "unitConversion.H"
#include "surfaceInterpolate.H"
#include "fvcDiv.H"
#include "fvcGrad.H"
#include "fvcSnGrad.H"


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// Correction for the boundary condition on the unit normal nHat on
// walls to produce the correct contact angle.

// The dynamic contact angle is calculated from the component of the
// velocity on the direction of the interface, parallel to the wall.

void Foam::interfaceProperties::correctContactAngle
(
    surfaceVectorField::Boundary& nHatb,
    const surfaceVectorField::Boundary& gradAlphaf
)
{
    const fvMesh& mesh = alpha1_.mesh();
    volScalarField::Boundary& a1bf = alpha1_.boundaryFieldRef();
    volScalarField::Boundary& a2bf = alpha2_.boundaryFieldRef();

    const fvBoundaryMesh& boundary = mesh.boundary();

    forAll(boundary, patchi)
    {
        if (isA<alphaContactAngleFvPatchScalarField>(a1bf[patchi]))
        {
            alphaContactAngleFvPatchScalarField& a1cap =
                refCast<alphaContactAngleFvPatchScalarField>
                (
                    a1bf[patchi]
                );

            fvsPatchVectorField& nHatp = nHatb[patchi];
            const scalarField theta
            (
                degToRad(a1cap.theta(U_.boundaryField()[patchi], nHatp))
            );

            const vectorField nf
            (
                boundary[patchi].nf()
            );

            // Reset nHatp to correspond to the contact angle

            const scalarField a12(nHatp & nf);
            const scalarField b1(cos(theta));

            scalarField b2(nHatp.size());
            forAll(b2, facei)
            {
                b2[facei] = cos(acos(a12[facei]) - theta[facei]);
            }

            const scalarField det(1.0 - a12*a12);

            scalarField a((b1 - a12*b2)/det);
            scalarField b((b2 - a12*b1)/det);

            nHatp = a*nf + b*nHatp;
            nHatp /= (mag(nHatp) + deltaN_.value());

            a1cap.gradient() = (nf & nHatp)*mag(gradAlphaf[patchi]);
            a1cap.evaluate();
            a2bf[patchi] = 1 - a1cap;
        }
    }
}


void Foam::interfaceProperties::calculateK()
{
    const fvMesh& mesh = alpha1_.mesh();
    const surfaceVectorField& Sf = mesh.Sf();

    // Cell gradient of alpha
    const volVectorField gradAlpha(fvc::grad(alpha1_, "nHat"));

    // Interpolated face-gradient of alpha
    surfaceVectorField gradAlphaf(fvc::interpolate(gradAlpha));

    // gradAlphaf -=
    //    (mesh.Sf()/mesh.magSf())
    //   *(fvc::snGrad(alpha1_) - (mesh.Sf() & gradAlphaf)/mesh.magSf());

    // Face unit interface normal
    surfaceVectorField nHatfv(gradAlphaf/(mag(gradAlphaf) + deltaN_));
    // surfaceVectorField nHatfv
    // (
    //     (gradAlphaf + deltaN_*vector(0, 0, 1)
    //    *sign(gradAlphaf.component(vector::Z)))/(mag(gradAlphaf) + deltaN_)
    // );
    correctContactAngle(nHatfv.boundaryFieldRef(), gradAlphaf.boundaryField());

    // Face unit interface normal flux
    nHatf_ = nHatfv & Sf;

    // Simple expression for curvature
    K_ = -fvc::div(nHatf_);

    // Complex expression for curvature.
    // Correction is formally zero but numerically non-zero.
    /*
    volVectorField nHat(gradAlpha/(mag(gradAlpha) + deltaN_));
    forAll(nHat.boundaryField(), patchi)
    {
        nHat.boundaryField()[patchi] = nHatfv.boundaryField()[patchi];
    }

    K_ = -fvc::div(nHatf_) + (nHat & fvc::grad(nHatfv) & nHat);
    */
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::interfaceProperties::interfaceProperties
(
    const IOdictionary& dict,
    volScalarField& alpha1,
    volScalarField& alpha2,
    const volVectorField& U
)
:
    phasePropertiesDict_(dict),

    alpha1_(alpha1),
    alpha2_(alpha2),
    U_(U),

    sigmaPtr_(surfaceTensionModel::New(dict, alpha1.mesh())),

    deltaN_
    (
        "deltaN",
        1e-8/pow(average(alpha1.mesh().V()), 1.0/3.0)
    ),

    nHatf_
    (
        IOobject
        (
            "nHatf",
            alpha1_.time().name(),
            alpha1_.mesh()
        ),
        alpha1_.mesh(),
        dimensionedScalar(dimArea, 0)
    ),

    K_
    (
        IOobject
        (
            "interfaceProperties:K",
            alpha1_.time().name(),
            alpha1_.mesh()
        ),
        alpha1_.mesh(),
        dimensionedScalar(dimless/dimLength, 0)
    )
{
    calculateK();
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::volVectorField> Foam::interfaceProperties::n() const
{
    const volVectorField gradAlpha(fvc::grad(alpha1_));

    return volVectorField::New("n", gradAlpha/(mag(gradAlpha) + deltaN_));
}


Foam::tmp<Foam::volScalarField>
Foam::interfaceProperties::sigmaK() const
{
    return sigmaPtr_->sigma()*K_;
}


Foam::tmp<Foam::surfaceScalarField>
Foam::interfaceProperties::surfaceTensionForce() const
{
    return fvc::interpolate(sigmaK())*fvc::snGrad(alpha1_);
}


Foam::tmp<Foam::volScalarField>
Foam::interfaceProperties::nearInterface() const
{
    return pos0(alpha1_ - 0.01)*pos0(0.99 - alpha1_);
}


Foam::tmp<Foam::volScalarField::Internal>
Foam::interfaceProperties::fraction() const
{
    const fvMesh& mesh = alpha1_.mesh();

    tmp<volScalarField::Internal> tA
    (
        volScalarField::Internal::New
        (
            "A",
            mesh,
            dimensionedScalar(dimless, 0)
        )
    );
    volScalarField::Internal& A = tA.ref();

    const surfaceVectorField& Sf = mesh.Sf();
    const labelUList& own = mesh.owner();
    const labelUList& nei = mesh.neighbour();

    const surfaceScalarField alphaf(fvc::interpolate(alpha1_));
    const volVectorField::Internal n(this->n());

    const scalarField& ialpha = alpha1_;
    const scalarField& ialphaf = alphaf;
    scalarField sumnSf(mesh.nCells(), 0);

    forAll(own, facei)
    {
        {
            const scalar nSf(mag(n[own[facei]] & Sf[facei]));
            A[own[facei]] += nSf*(ialphaf[facei] - ialpha[own[facei]]);
            sumnSf[own[facei]] += nSf;
        }
        {
            const scalar nSf(mag(n[nei[facei]] & Sf[facei]));
            A[nei[facei]] += nSf*(ialphaf[facei] - ialpha[nei[facei]]);
            sumnSf[nei[facei]] += nSf;
        }
    }

    forAll(mesh.boundary(), patchi)
    {
        const labelUList& own = mesh.boundary()[patchi].faceCells();
        const fvsPatchScalarField& palphaf = alphaf.boundaryField()[patchi];

        forAll(mesh.boundary()[patchi], facei)
        {
            const scalar nSf(mag(n[own[facei]] & Sf[facei]));
            A[own[facei]] += nSf*(palphaf[facei] - ialpha[own[facei]]);
            sumnSf[own[facei]] += nSf;
        }
    }

    scalarField& a = A.field();
    forAll(a, i)
    {
        if (sumnSf[i] > small)
        {
            a[i] = 2*mag(a[i])/sumnSf[i];
        }
        else
        {
            a[i] = 0;
        }
    }

    return tA;
}


void Foam::interfaceProperties::correct()
{
    calculateK();
}


bool Foam::interfaceProperties::read()
{
    sigmaPtr_->readDict(phasePropertiesDict_);

    return true;
}


// ************************************************************************* //
