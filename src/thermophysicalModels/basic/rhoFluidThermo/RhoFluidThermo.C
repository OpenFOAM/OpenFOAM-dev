/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2023 OpenFOAM Foundation
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

#include "RhoFluidThermo.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class BaseThermo>
void Foam::RhoFluidThermo<BaseThermo>::calculate()
{
    const scalarField& hCells = this->he();
    const scalarField& pCells = this->p_;

    scalarField& TCells = this->T_.primitiveFieldRef();
    scalarField& CpCells = this->Cp_.primitiveFieldRef();
    scalarField& CvCells = this->Cv_.primitiveFieldRef();
    scalarField& psiCells = this->psi_.primitiveFieldRef();
    scalarField& rhoCells = this->rho_.primitiveFieldRef();
    scalarField& muCells = this->mu_.primitiveFieldRef();
    scalarField& kappaCells = this->kappa_.primitiveFieldRef();

    auto Yslicer = this->Yslicer();

    forAll(TCells, celli)
    {
        auto composition = this->cellComposition(Yslicer, celli);

        const typename BaseThermo::mixtureType::thermoMixtureType&
            thermoMixture = this->thermoMixture(composition);

        const typename BaseThermo::mixtureType::transportMixtureType&
            transportMixture =
            this->transportMixture(composition, thermoMixture);

        TCells[celli] = thermoMixture.The
        (
            hCells[celli],
            pCells[celli],
            TCells[celli]
        );

        CpCells[celli] = thermoMixture.Cp(pCells[celli], TCells[celli]);
        CvCells[celli] = thermoMixture.Cv(pCells[celli], TCells[celli]);
        psiCells[celli] = thermoMixture.psi(pCells[celli], TCells[celli]);
        rhoCells[celli] = thermoMixture.rho(pCells[celli], TCells[celli]);

        muCells[celli] = transportMixture.mu(pCells[celli], TCells[celli]);
        kappaCells[celli] =
            transportMixture.kappa(pCells[celli], TCells[celli]);
    }

    volScalarField::Boundary& pBf =
        this->p_.boundaryFieldRef();

    volScalarField::Boundary& TBf =
        this->T_.boundaryFieldRef();

    volScalarField::Boundary& CpBf =
        this->Cp_.boundaryFieldRef();

    volScalarField::Boundary& CvBf =
        this->Cv_.boundaryFieldRef();

    volScalarField::Boundary& psiBf =
        this->psi_.boundaryFieldRef();

    volScalarField::Boundary& rhoBf =
        this->rho_.boundaryFieldRef();

    volScalarField::Boundary& heBf =
        this->he().boundaryFieldRef();

    volScalarField::Boundary& muBf =
        this->mu_.boundaryFieldRef();

    volScalarField::Boundary& kappaBf =
        this->kappa_.boundaryFieldRef();

    forAll(this->T_.boundaryField(), patchi)
    {
        fvPatchScalarField& pp = pBf[patchi];
        fvPatchScalarField& pT = TBf[patchi];
        fvPatchScalarField& pCp = CpBf[patchi];
        fvPatchScalarField& pCv = CvBf[patchi];
        fvPatchScalarField& ppsi = psiBf[patchi];
        fvPatchScalarField& prho = rhoBf[patchi];
        fvPatchScalarField& phe = heBf[patchi];
        fvPatchScalarField& pmu = muBf[patchi];
        fvPatchScalarField& pkappa = kappaBf[patchi];

        if (pT.fixesValue())
        {
            forAll(pT, facei)
            {
                auto composition =
                    this->patchFaceComposition(Yslicer, patchi, facei);

                const typename BaseThermo::mixtureType::thermoMixtureType&
                    thermoMixture = this->thermoMixture(composition);

                const typename BaseThermo::mixtureType::transportMixtureType&
                    transportMixture =
                    this->transportMixture(composition, thermoMixture);

                phe[facei] = thermoMixture.he(pp[facei], pT[facei]);

                pCp[facei] = thermoMixture.Cp(pp[facei], pT[facei]);
                pCv[facei] = thermoMixture.Cv(pp[facei], pT[facei]);
                ppsi[facei] = thermoMixture.psi(pp[facei], pT[facei]);
                prho[facei] = thermoMixture.rho(pp[facei], pT[facei]);

                pmu[facei] = transportMixture.mu(pp[facei], pT[facei]);
                pkappa[facei] = transportMixture.kappa(pp[facei], pT[facei]);
            }
        }
        else
        {
            forAll(pT, facei)
            {
                auto composition =
                    this->patchFaceComposition(Yslicer, patchi, facei);

                const typename BaseThermo::mixtureType::thermoMixtureType&
                    thermoMixture = this->thermoMixture(composition);

                const typename BaseThermo::mixtureType::transportMixtureType&
                    transportMixture =
                    this->transportMixture(composition, thermoMixture);

                pT[facei] = thermoMixture.The(phe[facei], pp[facei], pT[facei]);

                pCp[facei] = thermoMixture.Cp(pp[facei], pT[facei]);
                pCv[facei] = thermoMixture.Cv(pp[facei], pT[facei]);
                ppsi[facei] = thermoMixture.psi(pp[facei], pT[facei]);
                prho[facei] = thermoMixture.rho(pp[facei], pT[facei]);

                pmu[facei] = transportMixture.mu(pp[facei], pT[facei]);
                pkappa[facei] = transportMixture.kappa(pp[facei], pT[facei]);
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BaseThermo>
Foam::RhoFluidThermo<BaseThermo>::RhoFluidThermo
(
    const fvMesh& mesh,
    const word& phaseName
)
:
    BaseThermo(mesh, phaseName)
{
    calculate();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class BaseThermo>
Foam::RhoFluidThermo<BaseThermo>::~RhoFluidThermo()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BaseThermo>
void Foam::RhoFluidThermo<BaseThermo>::correct()
{
    if (BaseThermo::debug)
    {
        InfoInFunction << endl;
    }

    calculate();

    if (BaseThermo::debug)
    {
        Info<< "    Finished" << endl;
    }
}


// ************************************************************************* //
