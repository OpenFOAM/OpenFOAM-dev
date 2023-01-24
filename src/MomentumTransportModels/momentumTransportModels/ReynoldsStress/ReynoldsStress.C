/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2015-2023 OpenFOAM Foundation
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

#include "ReynoldsStress.H"
#include "fvcSnGrad.H"
#include "wallFvPatch.H"

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class BasicMomentumTransportModel>
void Foam::ReynoldsStress<BasicMomentumTransportModel>::boundNormalStress
(
    volSymmTensorField& R
) const
{
    scalar kMin = this->kMin_.value();

    R.max
    (
        dimensionedSymmTensor
        (
            "zero",
            R.dimensions(),
            symmTensor
            (
                kMin, -great, -great,
                kMin, -great,
                kMin
            )
        )
    );
}


template<class BasicMomentumTransportModel>
void Foam::ReynoldsStress<BasicMomentumTransportModel>::correctWallShearStress
(
    volSymmTensorField& R
) const
{
    const fvPatchList& patches = this->mesh_.boundary();

    volSymmTensorField::Boundary& RBf = R.boundaryFieldRef();

    forAll(patches, patchi)
    {
        const fvPatch& curPatch = patches[patchi];

        if (isA<wallFvPatch>(curPatch))
        {
            symmTensorField& Rw = RBf[patchi];

            const scalarField& nutw = this->nut_.boundaryField()[patchi];

            const vectorField snGradUw
            (
                this->U_.boundaryField()[patchi].snGrad()
            );

            const vectorField& Sf = this->mesh_.Sf().boundaryField()[patchi];

            const scalarField& magSf =
                this->mesh_.magSf().boundaryField()[patchi];

            forAll(curPatch, facei)
            {
                // Calculate near-wall velocity gradient
                const tensor gradUw
                    = (Sf[facei]/magSf[facei])*snGradUw[facei];

                // Set the wall Reynolds-stress to the near-wall shear-stress
                // Note: the spherical part of the normal stress is included in
                // the pressure
                Rw[facei] = -nutw[facei]*dev(twoSymm(gradUw));
            }
        }
    }
}


template<class BasicMomentumTransportModel>
Foam::tmp<Foam::fvSymmTensorMatrix>
Foam::ReynoldsStress<BasicMomentumTransportModel>::RSource() const
{
    return tmp<fvSymmTensorMatrix>
    (
        new fvSymmTensorMatrix
        (
            this->R_,
            dimVolume*this->rho_.dimensions()*this->R_.dimensions()/dimTime
        )
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicMomentumTransportModel>
Foam::ReynoldsStress<BasicMomentumTransportModel>::ReynoldsStress
(
    const word& modelName,
    const alphaField& alpha,
    const rhoField& rho,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const viscosity& viscosity
)
:
    BasicMomentumTransportModel
    (
        modelName,
        alpha,
        rho,
        U,
        alphaRhoPhi,
        phi,
        viscosity
    ),

    couplingFactor_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "couplingFactor",
            this->coeffDict_,
            1
        )
    ),

    R_
    (
        IOobject
        (
            this->groupName("R"),
            this->runTime_.name(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    ),

    nut_
    (
        IOobject
        (
            this->groupName("nut"),
            this->runTime_.name(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    )
{
    if (couplingFactor_.value() < 0 || couplingFactor_.value() > 1)
    {
        FatalErrorInFunction
            << "couplingFactor = " << couplingFactor_
            << " is not in range 0 - 1" << nl
            << exit(FatalError);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicMomentumTransportModel>
bool Foam::ReynoldsStress<BasicMomentumTransportModel>::read()
{
    return BasicMomentumTransportModel::read();
}


template<class BasicMomentumTransportModel>
Foam::tmp<Foam::volSymmTensorField>
Foam::ReynoldsStress<BasicMomentumTransportModel>::sigma() const
{
    return R_;
}


template<class BasicMomentumTransportModel>
Foam::tmp<Foam::volScalarField>
Foam::ReynoldsStress<BasicMomentumTransportModel>::k() const
{
    tmp<Foam::volScalarField> tk(0.5*tr(R_));
    tk.ref().rename("k");
    return tk;
}


template<class BasicMomentumTransportModel>
Foam::tmp<Foam::volSymmTensorField>
Foam::ReynoldsStress<BasicMomentumTransportModel>::devTau() const
{
    return volSymmTensorField::New
    (
        this->groupName("devTau"),
        this->alpha_*this->rho_*R_
      - (this->alpha_*this->rho_*this->nu())
       *dev(twoSymm(fvc::grad(this->U_)))
    );
}


template<class BasicMomentumTransportModel>
template<class RhoFieldType>
Foam::tmp<Foam::fvVectorMatrix>
Foam::ReynoldsStress<BasicMomentumTransportModel>::DivDevRhoReff
(
    const RhoFieldType& rho,
    volVectorField& U
) const
{
    tmp<volTensorField> tgradU = fvc::grad(U);
    const volTensorField& gradU = tgradU();
    const surfaceTensorField gradUf(fvc::interpolate(gradU));

    // Interpolate Reynolds stress to the faces
    // with either a stress or velocity coupling correction
    const surfaceVectorField Refff
    (
        (this->mesh().Sf() & fvc::interpolate(R_))

        // Stress coupling
      + couplingFactor_
       *(this->mesh().Sf() & fvc::interpolate(this->nut()*gradU))

        // or velocity gradient coupling
   // + couplingFactor_
   //  *fvc::interpolate(this->nut())*(this->mesh().Sf() & gradUf)

      - fvc::interpolate(couplingFactor_*this->nut() + this->nu())
       *this->mesh().magSf()*fvc::snGrad(U)
      - fvc::interpolate(this->nu())*(this->mesh().Sf() & dev2(gradUf.T()))
    );

    return
    (
        fvc::div(fvc::interpolate(this->alpha_*rho)*Refff)
      - correction(fvm::laplacian(this->alpha_*rho*this->nuEff(), U))
    );
}


template<class BasicMomentumTransportModel>
Foam::tmp<Foam::fvVectorMatrix>
Foam::ReynoldsStress<BasicMomentumTransportModel>::divDevTau
(
    volVectorField& U
) const
{
    return DivDevRhoReff(this->rho_, U);
}


template<class BasicMomentumTransportModel>
Foam::tmp<Foam::fvVectorMatrix>
Foam::ReynoldsStress<BasicMomentumTransportModel>::divDevTau
(
    const volScalarField& rho,
    volVectorField& U
) const
{
    return DivDevRhoReff(rho, U);
}


template<class BasicMomentumTransportModel>
void Foam::ReynoldsStress<BasicMomentumTransportModel>::validate()
{
    correctNut();
}


template<class BasicMomentumTransportModel>
void Foam::ReynoldsStress<BasicMomentumTransportModel>::correct()
{
    BasicMomentumTransportModel::correct();
}


// ************************************************************************* //
