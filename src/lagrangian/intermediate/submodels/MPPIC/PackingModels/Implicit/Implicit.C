/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013-2014 OpenFOAM Foundation
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

#include "Implicit.H"
#include "fixedValueFvsPatchField.H"
#include "fvmDdt.H"
#include "fvmDiv.H"
#include "fvmLaplacian.H"
#include "fvcReconstruct.H"
#include "volPointInterpolation.H"                                              

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::PackingModels::Implicit<CloudType>::Implicit
(
    const dictionary& dict,
    CloudType& owner
)
:
    PackingModel<CloudType>(dict, owner, typeName),
    alpha_
    (
        this->owner().name() + ":alpha",
        this->owner().theta()
    ),
    phiCorrect_(NULL),
    uCorrect_(NULL),
    applyGravity_(this->coeffDict().lookup("applyGravity")),
    alphaMin_(readScalar(this->coeffDict().lookup("alphaMin"))),
    rhoMin_(readScalar(this->coeffDict().lookup("rhoMin")))
{
    alpha_.oldTime();
}


template<class CloudType>
Foam::PackingModels::Implicit<CloudType>::Implicit
(
    const Implicit<CloudType>& cm
)
:
    PackingModel<CloudType>(cm),
    alpha_(cm.alpha_),
    phiCorrect_(cm.phiCorrect_()),
    uCorrect_(cm.uCorrect_()),
    applyGravity_(cm.applyGravity_),
    alphaMin_(cm.alphaMin_),
    rhoMin_(cm.rhoMin_)
{
    alpha_.oldTime();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::PackingModels::Implicit<CloudType>::~Implicit()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
void Foam::PackingModels::Implicit<CloudType>::cacheFields(const bool store)
{
    PackingModel<CloudType>::cacheFields(store);

    if (store)
    { 
        const fvMesh& mesh = this->owner().mesh();
        const dimensionedScalar deltaT = this->owner().db().time().deltaT();
        const word& cloudName = this->owner().name();

        const dimensionedVector& g = this->owner().g();
        const volScalarField& rhoc = this->owner().rho();

        const AveragingMethod<scalar>& rhoAverage =
            mesh.lookupObject<AveragingMethod<scalar> >
            (
                cloudName + ":rhoAverage"
            );
        const AveragingMethod<scalar>& uSqrAverage =
            mesh.lookupObject<AveragingMethod<scalar> >
            (
                cloudName + ":uSqrAverage"
            );


        // Property fields
        // ~~~~~~~~~~~~~~~

        // volume fraction field
        alpha_ = max(this->owner().theta(), alphaMin_);
        alpha_.correctBoundaryConditions();

        // average density
        volScalarField rho
        (
            IOobject
            (
                cloudName + ":rho",
                this->owner().db().time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar("zero", dimDensity, 0),
            zeroGradientFvPatchField<scalar>::typeName
        );
        rho.internalField() = max(rhoAverage.internalField(), rhoMin_);
        rho.correctBoundaryConditions();

        //Info << "       x: " << mesh.C().internalField().component(2) << endl;
        //Info << "   alpha: " << alpha_.internalField() << endl;
        //Info << "alphaOld: " << alpha_.oldTime().internalField() << endl;
        //Info << "     rho: " << rho.internalField() << endl;
        //Info << endl;

        
        // Stress field
        // ~~~~~~~~~~~~

        // stress derivative wrt volume fraction
        volScalarField tauPrime
        (
            IOobject
            (
                cloudName + ":tauPrime",
                this->owner().db().time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar("zero", dimPressure, 0),
            zeroGradientFvPatchField<scalar>::typeName
        );

        tauPrime.internalField() =
            this->particleStressModel_->dTaudTheta
            (
                alpha_.internalField(),
                rho.internalField(),
                uSqrAverage.internalField()
            )();

        tauPrime.correctBoundaryConditions();


        // Implicit solution for the volume fraction
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        surfaceScalarField
            tauPrimeByRhoAf
            (
                "tauPrimeByRhoAf",
                fvc::interpolate(deltaT*tauPrime/rho)
            );

        fvScalarMatrix alphaEqn
        (
            fvm::ddt(alpha_)
          - fvc::ddt(alpha_)
          - fvm::laplacian(tauPrimeByRhoAf, alpha_)
        );

        if (applyGravity_)
        {
            surfaceScalarField
                phiGByA
                (
                    "phiGByA",
                    deltaT*(g & mesh.Sf())*fvc::interpolate(1.0 - rhoc/rho)
                );

            alphaEqn += fvm::div(phiGByA, alpha_);
        }

        alphaEqn.solve();


        // Generate correction fields
        // ~~~~~~~~~~~~~~~~~

        // correction volumetric flux
        phiCorrect_ = tmp<surfaceScalarField>
        (
            new surfaceScalarField
            (
                cloudName + ":phiCorrect",
                alphaEqn.flux()/fvc::interpolate(alpha_)
            )
        );

        // correction velocity
        uCorrect_ = tmp<volVectorField>
        (
            new volVectorField
            (
                cloudName + ":uCorrect",
                fvc::reconstruct(phiCorrect_())
            )

        );
        uCorrect_->correctBoundaryConditions();

        //Info << endl;
        //Info << "     alpha: " << alpha_.internalField() << endl;
        //Info << "phiCorrect: " << phiCorrect_->internalField() << endl;
        //Info << "  uCorrect: " << uCorrect_->internalField() << endl;
        //Info << endl;
    }
    else
    {
        alpha_.oldTime();
        phiCorrect_.clear();
        uCorrect_.clear();
    }
}


template<class CloudType>
Foam::vector Foam::PackingModels::Implicit<CloudType>::velocityCorrection
(
    typename CloudType::parcelType& p,
    const scalar deltaT
) const
{
    const fvMesh& mesh = this->owner().mesh();

    // containing tetrahedron and parcel coordinates within
    const label cellI = p.cell();
    const label faceI = p.tetFace();
    const tetIndices tetIs(cellI, faceI, p.tetPt(), mesh);
    List<scalar> tetCoordinates(4);
    tetIs.tet(mesh).barycentric(p.position(), tetCoordinates);

    // cell velocity
    const vector U = uCorrect_()[cellI];

    // face geometry
    vector nHat = mesh.faces()[faceI].normal(mesh.points());
    const scalar nMag = mag(nHat);
    nHat /= nMag;

    // get face flux
    scalar phi;
    const label patchI = mesh.boundaryMesh().whichPatch(faceI);
    if (patchI == -1)
    {
        phi = phiCorrect_()[faceI];
    }
    else
    {
        phi =
            phiCorrect_().boundaryField()[patchI]
            [
                mesh.boundaryMesh()[patchI].whichFace(faceI)
            ];
    }

    // interpolant equal to 1 at the cell centre and 0 at the face
    const scalar t = tetCoordinates[0];

    // the normal component of the velocity correction is interpolated linearly
    // the tangential component is equal to that at the cell centre
    return U + (1.0 - t)*nHat*(phi/nMag - (U & nHat));
}


// ************************************************************************* //
