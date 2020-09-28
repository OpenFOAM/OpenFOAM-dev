/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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

#include "heSolidThermo.H"
#include "fvmLaplacian.H"
#include "fvcLaplacian.H"
#include "coordinateSystem.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class BasicSolidThermo, class MixtureType>
void Foam::heSolidThermo<BasicSolidThermo, MixtureType>::calculate()
{
    const scalarField& hCells = this->he_;
    const auto& pCells = this->p_;

    scalarField& TCells = this->T_.primitiveFieldRef();
    scalarField& rhoCells = this->rho_.primitiveFieldRef();
    scalarField& alphaCells = this->alpha_.primitiveFieldRef();

    forAll(TCells, celli)
    {
        const typename MixtureType::thermoMixtureType& thermoMixture =
            this->cellThermoMixture(celli);

        const typename MixtureType::transportMixtureType& transportMixture =
            this->cellTransportMixture(celli, thermoMixture);

        rhoCells[celli] = thermoMixture.rho(pCells[celli], TCells[celli]);

        TCells[celli] = thermoMixture.THE
        (
            hCells[celli],
            pCells[celli],
            TCells[celli]
        );

        alphaCells[celli] =
            transportMixture.kappa(pCells[celli], TCells[celli])
           /thermoMixture.Cpv(pCells[celli], TCells[celli]);
    }

    const auto& pBf = this->p_.boundaryField();

    volScalarField::Boundary& TBf =
        this->T_.boundaryFieldRef();

    volScalarField::Boundary& rhoBf =
        this->rho_.boundaryFieldRef();

    volScalarField::Boundary& heBf =
        this->he().boundaryFieldRef();

    volScalarField::Boundary& alphaBf =
        this->alpha_.boundaryFieldRef();

    forAll(this->T_.boundaryField(), patchi)
    {
        const auto& pp = pBf[patchi];
        fvPatchScalarField& pT = TBf[patchi];
        fvPatchScalarField& prho = rhoBf[patchi];
        fvPatchScalarField& phe = heBf[patchi];
        fvPatchScalarField& palpha = alphaBf[patchi];

        if (pT.fixesValue())
        {
            forAll(pT, facei)
            {
                const typename MixtureType::thermoMixtureType&
                    thermoMixture = this->patchFaceThermoMixture(patchi, facei);

                const typename MixtureType::transportMixtureType&
                    transportMixture =
                    this->patchFaceTransportMixture
                    (patchi, facei, thermoMixture);

                phe[facei] = thermoMixture.HE(pp[facei], pT[facei]);
                prho[facei] = thermoMixture.rho(pp[facei], pT[facei]);

                palpha[facei] =
                    transportMixture.kappa(pp[facei], pT[facei])
                   /thermoMixture.Cpv(pp[facei], pT[facei]);
            }
        }
        else
        {
            forAll(pT, facei)
            {
                const typename MixtureType::thermoMixtureType&
                    thermoMixture = this->patchFaceThermoMixture(patchi, facei);

                const typename MixtureType::transportMixtureType&
                    transportMixture =
                    this->patchFaceTransportMixture
                    (patchi, facei, thermoMixture);

                pT[facei] = thermoMixture.THE(phe[facei], pp[facei] ,pT[facei]);
                prho[facei] = thermoMixture.rho(pp[facei], pT[facei]);

                palpha[facei] =
                    transportMixture.kappa(pp[facei], pT[facei])
                   /thermoMixture.Cpv(pp[facei], pT[facei]);
            }
        }
    }

    this->alpha_.correctBoundaryConditions();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicSolidThermo, class MixtureType>
Foam::heSolidThermo<BasicSolidThermo, MixtureType>::
heSolidThermo
(
    const fvMesh& mesh,
    const word& phaseName
)
:
    heThermo<BasicSolidThermo, MixtureType>(mesh, phaseName)
{
    calculate();
}


template<class BasicSolidThermo, class MixtureType>
Foam::heSolidThermo<BasicSolidThermo, MixtureType>::
heSolidThermo
(
    const fvMesh& mesh,
    const dictionary& dict,
    const word& phaseName
)
:
    heThermo<BasicSolidThermo, MixtureType>(mesh, dict, phaseName)
{
    calculate();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class BasicSolidThermo, class MixtureType>
Foam::heSolidThermo<BasicSolidThermo, MixtureType>::~heSolidThermo()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicSolidThermo, class MixtureType>
void Foam::heSolidThermo<BasicSolidThermo, MixtureType>::correct()
{
    if (debug)
    {
        InfoInFunction << endl;
    }

    calculate();

    if (debug)
    {
        Info<< "    Finished" << endl;
    }
}


template<class BasicSolidThermo, class MixtureType>
Foam::tmp<Foam::volVectorField>
Foam::heSolidThermo<BasicSolidThermo, MixtureType>::Kappa() const
{
    const fvMesh& mesh = this->T_.mesh();

    tmp<volVectorField> tKappa
    (
        volVectorField::New
        (
            "Kappa",
            mesh,
            dimEnergy/dimTime/dimLength/dimTemperature
        )
    );

    const auto& pCells = this->p_;
    const scalarField& TCells = this->T_;

    volVectorField& Kappa = tKappa.ref();
    vectorField& KappaCells = Kappa.primitiveFieldRef();

    forAll(KappaCells, celli)
    {
        Kappa[celli] =
            this->cellTransportMixture
            (
                celli
            ).Kappa(pCells[celli], TCells[celli]);
    }

    volVectorField::Boundary& KappaBf = Kappa.boundaryFieldRef();

    forAll(KappaBf, patchi)
    {
        const auto& pp = this->p_.boundaryField()[patchi];
        const scalarField& pT = this->T_.boundaryField()[patchi];

        vectorField& Kappap = KappaBf[patchi];

        forAll(Kappap, facei)
        {
            Kappap[facei] =
                this->patchFaceTransportMixture
                (
                    patchi,
                    facei
                ).Kappa(pp[facei], pT[facei]);
        }
    }

    return tKappa;
}


template<class BasicSolidThermo, class MixtureType>
Foam::tmp<Foam::vectorField>
Foam::heSolidThermo<BasicSolidThermo, MixtureType>::Kappa
(
    const label patchi
) const
{
    const auto& pp = this->p_.boundaryField()[patchi];
    const scalarField& Tp = this->T_.boundaryField()[patchi];

    tmp<vectorField> tKappa(new vectorField(Tp.size()));
    vectorField& Kappap = tKappa.ref();

    forAll(Tp, facei)
    {
        Kappap[facei] =
            this->patchFaceTransportMixture
            (
                patchi,
                facei
            ).Kappa(pp[patchi], Tp[facei]);
    }

    return tKappa;
}


template<class BasicSolidThermo, class MixtureType>
Foam::tmp<Foam::volSymmTensorField>
Foam::heSolidThermo<BasicSolidThermo, MixtureType>::KappaLocal() const
{
    const fvMesh& mesh = this->T_.mesh();

    const coordinateSystem coordinates
    (
        coordinateSystem::New(mesh, this->properties())
    );

    const tmp<volVectorField> tKappa(Kappa());
    const volVectorField& Kappa = tKappa();

    tmp<volSymmTensorField> tKappaLocal
    (
        volSymmTensorField::New
        (
            "KappaLocal",
            mesh,
            dimensionedSymmTensor(Kappa.dimensions(), Zero)
        )
    );
    volSymmTensorField& KappaLocal = tKappaLocal.ref();

    KappaLocal.primitiveFieldRef() =
        coordinates.R(mesh.C()).transformVector(Kappa);

    forAll(KappaLocal.boundaryField(), patchi)
    {
        KappaLocal.boundaryFieldRef()[patchi] =
            coordinates.R(mesh.boundary()[patchi].Cf())
           .transformVector(Kappa.boundaryField()[patchi]);
    }

    return tKappaLocal;
}


template<class BasicSolidThermo, class MixtureType>
Foam::tmp<Foam::symmTensorField>
Foam::heSolidThermo<BasicSolidThermo, MixtureType>::KappaLocal
(
    const label patchi
) const
{
    const fvMesh& mesh = this->T_.mesh();

    const coordinateSystem coordinates
    (
        coordinateSystem::New(mesh, this->properties())
    );

    return
        coordinates.R(mesh.boundary()[patchi].Cf())
       .transformVector(Kappa(patchi));
}


template<class BasicSolidThermo, class MixtureType>
Foam::tmp<Foam::surfaceScalarField>
Foam::heSolidThermo<BasicSolidThermo, MixtureType>::q() const
{
    const fvMesh& mesh = this->T_.mesh();
    mesh.setFluxRequired(this->T_.name());

    return
      - (
            isotropic()
          ? fvm::laplacian(this->kappa(), this->T_)().flux()
          : fvm::laplacian(KappaLocal(), this->T_)().flux()
        );
}


template<class BasicSolidThermo, class MixtureType>
Foam::tmp<Foam::fvScalarMatrix>
Foam::heSolidThermo<BasicSolidThermo, MixtureType>::divq
(
    volScalarField& e
) const
{
    return
      - (
            isotropic()
          ?   fvc::laplacian(this->kappa(), this->T_)
            + correction(fvm::laplacian(this->alpha(), e))
          :   fvc::laplacian(KappaLocal(), this->T_)
            + correction
              (
                  fvm::laplacian
                  (
                      KappaLocal()/this->Cv(),
                      e,
                      "laplacian(" + this->alpha().name() + ",e)"
                  )
              )
        );
}


// ************************************************************************* //
