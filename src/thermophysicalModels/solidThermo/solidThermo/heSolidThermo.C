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
    scalarField& CpCells = this->Cp_.primitiveFieldRef();
    scalarField& CvCells = this->Cv_.primitiveFieldRef();
    scalarField& rhoCells = this->rho_.primitiveFieldRef();
    scalarField& kappaCells = this->kappa_.primitiveFieldRef();

    forAll(TCells, celli)
    {
        const typename MixtureType::thermoMixtureType& thermoMixture =
            this->cellThermoMixture(celli);

        const typename MixtureType::transportMixtureType& transportMixture =
            this->cellTransportMixture(celli, thermoMixture);

        TCells[celli] = thermoMixture.THE
        (
            hCells[celli],
            pCells[celli],
            TCells[celli]
        );

        CpCells[celli] = thermoMixture.Cp(pCells[celli], TCells[celli]);
        CvCells[celli] = thermoMixture.Cv(pCells[celli], TCells[celli]);
        rhoCells[celli] = thermoMixture.rho(pCells[celli], TCells[celli]);

        kappaCells[celli] =
            transportMixture.kappa(pCells[celli], TCells[celli]);
    }


    volScalarField::Boundary& heBf =
        this->he().boundaryFieldRef();

    const auto& pBf = this->p_.boundaryField();

    volScalarField::Boundary& TBf =
        this->T_.boundaryFieldRef();

    volScalarField::Boundary& CpBf =
        this->Cp_.boundaryFieldRef();

    volScalarField::Boundary& CvBf =
        this->Cv_.boundaryFieldRef();

    volScalarField::Boundary& rhoBf =
        this->rho_.boundaryFieldRef();

    volScalarField::Boundary& kappaBf =
        this->kappa_.boundaryFieldRef();

    forAll(this->T_.boundaryField(), patchi)
    {
        fvPatchScalarField& phe = heBf[patchi];
        const auto& pp = pBf[patchi];
        fvPatchScalarField& pT = TBf[patchi];
        fvPatchScalarField& pCp = CpBf[patchi];
        fvPatchScalarField& pCv = CvBf[patchi];
        fvPatchScalarField& prho = rhoBf[patchi];
        fvPatchScalarField& pkappa = kappaBf[patchi];

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
                pCp[facei] = thermoMixture.Cp(pp[facei], pT[facei]);
                pCv[facei] = thermoMixture.Cv(pp[facei], pT[facei]);

                pkappa[facei] = transportMixture.kappa(pp[facei], pT[facei]);
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
                pCp[facei] = thermoMixture.Cp(pp[facei], pT[facei]);
                pCv[facei] = thermoMixture.Cv(pp[facei], pT[facei]);

                pkappa[facei] = transportMixture.kappa(pp[facei], pT[facei]);
            }
        }
    }

    this->kappa_.correctBoundaryConditions();
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
    const fvMesh& mesh = this->mesh();

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


// ************************************************************************* //
