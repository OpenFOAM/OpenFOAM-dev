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

#include "laminar.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "fvcGrad.H"
#include "fvcDiv.H"
#include "fvmLaplacian.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
Foam::laminar<BasicTurbulenceModel>::laminar
(
    const alphaField& alpha,
    const rhoField& rho,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const transportModel& transport,
    const word& propertiesName
)
:
    BasicTurbulenceModel
    (
        alpha,
        rho,
        U,
        alphaRhoPhi,
        phi,
        transport,
        propertiesName
    )
{}


// * * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
Foam::autoPtr<Foam::laminar<BasicTurbulenceModel> >
Foam::laminar<BasicTurbulenceModel>::New
(
    const alphaField& alpha,
    const rhoField& rho,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const transportModel& transport,
    const word& propertiesName
)
{
    return autoPtr<laminar>
    (
        new laminar
        (
            alpha,
            rho,
            U,
            alphaRhoPhi,
            phi,
            transport,
            propertiesName
        )
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
const Foam::dictionary&
Foam::laminar<BasicTurbulenceModel>::coeffDict() const
{
    return dictionary::null;
}


template<class BasicTurbulenceModel>
Foam::tmp<Foam::volScalarField>
Foam::laminar<BasicTurbulenceModel>::nut() const
{
    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                IOobject::groupName("nut", this->U_.group()),
                this->runTime_.timeName(),
                this->mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            this->mesh_,
            dimensionedScalar("nut", dimViscosity, 0.0)
        )
    );
}


template<class BasicTurbulenceModel>
Foam::tmp<Foam::scalarField>
Foam::laminar<BasicTurbulenceModel>::nut
(
    const label patchi
) const
{
    return tmp<scalarField>
    (
        new scalarField(this->mesh_.boundary()[patchi].size(), 0.0)
    );
}


template<class BasicTurbulenceModel>
Foam::tmp<Foam::volScalarField>
Foam::laminar<BasicTurbulenceModel>::nuEff() const
{
    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject::groupName("nuEff", this->U_.group()), this->nu()
        )
    );
}


template<class BasicTurbulenceModel>
Foam::tmp<Foam::scalarField>
Foam::laminar<BasicTurbulenceModel>::nuEff
(
    const label patchi
) const
{
    return this->nu(patchi);
}


template<class BasicTurbulenceModel>
Foam::tmp<Foam::volScalarField>
Foam::laminar<BasicTurbulenceModel>::k() const
{
    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                IOobject::groupName("k", this->U_.group()),
                this->runTime_.timeName(),
                this->mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            this->mesh_,
            dimensionedScalar("k", sqr(this->U_.dimensions()), 0.0)
        )
    );
}


template<class BasicTurbulenceModel>
Foam::tmp<Foam::volScalarField>
Foam::laminar<BasicTurbulenceModel>::epsilon() const
{
    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                IOobject::groupName("epsilon", this->U_.group()),
                this->runTime_.timeName(),
                this->mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            this->mesh_,
            dimensionedScalar
            (
                "epsilon", sqr(this->U_.dimensions())/dimTime, 0.0
            )
        )
    );
}


template<class BasicTurbulenceModel>
Foam::tmp<Foam::volSymmTensorField>
Foam::laminar<BasicTurbulenceModel>::R() const
{
    return tmp<volSymmTensorField>
    (
        new volSymmTensorField
        (
            IOobject
            (
                IOobject::groupName("R", this->U_.group()),
                this->runTime_.timeName(),
                this->mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            this->mesh_,
            dimensionedSymmTensor
            (
                "R", sqr(this->U_.dimensions()), symmTensor::zero
            )
        )
    );
}


template<class BasicTurbulenceModel>
Foam::tmp<Foam::volSymmTensorField>
Foam::laminar<BasicTurbulenceModel>::devRhoReff() const
{
    return tmp<volSymmTensorField>
    (
        new volSymmTensorField
        (
            IOobject
            (
                IOobject::groupName("devRhoReff", this->U_.group()),
                this->runTime_.timeName(),
                this->mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
           -(this->alpha_*this->rho_*nuEff())*dev(twoSymm(fvc::grad(this->U_)))
        )
    );
}


template<class BasicTurbulenceModel>
Foam::tmp<Foam::fvVectorMatrix>
Foam::laminar<BasicTurbulenceModel>::divDevRhoReff
(
    volVectorField& U
) const
{
    return
    (
      - fvm::laplacian(this->alpha_*this->rho_*nuEff(), U)
      - fvc::div(this->alpha_*this->rho_*nuEff()*dev2(T(fvc::grad(U))))
    );
}


template<class BasicTurbulenceModel>
void Foam::laminar<BasicTurbulenceModel>::correct()
{
    BasicTurbulenceModel::correct();
}


template<class BasicTurbulenceModel>
bool Foam::laminar<BasicTurbulenceModel>::read()
{
    return true;
}


// ************************************************************************* //
