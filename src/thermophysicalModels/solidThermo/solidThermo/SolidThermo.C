/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2024 OpenFOAM Foundation
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

#include "SolidThermo.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class BaseThermo>
void Foam::SolidThermo<BaseThermo>::calculate()
{
    const bool isotropic = this->isotropic();

    const scalarField& hCells = this->he_;
    const auto& pCells = this->p_;

    scalarField& TCells = this->T_.primitiveFieldRef();
    scalarField& CpCells = this->Cp_.primitiveFieldRef();
    scalarField& CvCells = this->Cv_.primitiveFieldRef();
    scalarField& rhoCells = this->rho_.primitiveFieldRef();
    scalarField& kappaCells = this->kappa_.primitiveFieldRef();
    vectorField& KappaCells = this->Kappa_.primitiveFieldRef();

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
        rhoCells[celli] = thermoMixture.rho(pCells[celli], TCells[celli]);

        if (isotropic)
        {
            kappaCells[celli] =
                transportMixture.kappa(pCells[celli], TCells[celli]);
        }
        else
        {
            KappaCells[celli] =
                transportMixture.Kappa(pCells[celli], TCells[celli]);
        }
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

    volVectorField::Boundary& KappaBf =
        this->Kappa_.boundaryFieldRef();

    forAll(this->T_.boundaryField(), patchi)
    {
        fvPatchScalarField& phe = heBf[patchi];
        const auto& pp = pBf[patchi];
        fvPatchScalarField& pT = TBf[patchi];
        fvPatchScalarField& pCp = CpBf[patchi];
        fvPatchScalarField& pCv = CvBf[patchi];
        fvPatchScalarField& prho = rhoBf[patchi];
        fvPatchScalarField& pkappa = kappaBf[patchi];
        fvPatchVectorField& pKappa = KappaBf[patchi];

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

                prho[facei] = thermoMixture.rho(pp[facei], pT[facei]);
                pCp[facei] = thermoMixture.Cp(pp[facei], pT[facei]);
                pCv[facei] = thermoMixture.Cv(pp[facei], pT[facei]);

                if (isotropic)
                {
                    pkappa[facei] =
                        transportMixture.kappa(pp[facei], pT[facei]);
                }
                else
                {
                    pKappa[facei] =
                        transportMixture.Kappa(pp[facei], pT[facei]);
                }
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

                pT[facei] = thermoMixture.The(phe[facei], pp[facei] ,pT[facei]);

                prho[facei] = thermoMixture.rho(pp[facei], pT[facei]);
                pCp[facei] = thermoMixture.Cp(pp[facei], pT[facei]);
                pCv[facei] = thermoMixture.Cv(pp[facei], pT[facei]);

                if (isotropic)
                {
                    pkappa[facei] =
                        transportMixture.kappa(pp[facei], pT[facei]);
                }
                else
                {
                    pKappa[facei] =
                        transportMixture.Kappa(pp[facei], pT[facei]);
                }
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BaseThermo>
Foam::SolidThermo<BaseThermo>::SolidThermo
(
    const fvMesh& mesh,
    const word& phaseName
)
:
    BaseThermo(mesh, phaseName),
    Kappa_
    (
        IOobject
        (
            BaseThermo::phasePropertyName("Kappa", phaseName),
            mesh.time().name(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedVector(dimThermalConductivity, Zero)
    )
{
    calculate();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class BaseThermo>
Foam::SolidThermo<BaseThermo>::~SolidThermo()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BaseThermo>
void Foam::SolidThermo<BaseThermo>::correct()
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


template<class BaseThermo>
const Foam::volVectorField&
Foam::SolidThermo<BaseThermo>::Kappa() const
{
    return Kappa_;
}


// ************************************************************************* //
