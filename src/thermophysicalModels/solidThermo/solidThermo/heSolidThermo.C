/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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
#include "volFields.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class BasicSolidThermo, class MixtureType>
void Foam::heSolidThermo<BasicSolidThermo, MixtureType>::calculate()
{
    scalarField& TCells = this->T_.primitiveFieldRef();

    const scalarField& hCells = this->he_;
    const scalarField& pCells = this->p_;
    scalarField& rhoCells = this->rho_.primitiveFieldRef();
    scalarField& alphaCells = this->alpha_.primitiveFieldRef();

    forAll(TCells, celli)
    {
        const typename MixtureType::thermoType& mixture_ =
            this->cellMixture(celli);

        const typename MixtureType::thermoType& volMixture_ =
            this->cellVolMixture(pCells[celli], TCells[celli], celli);

        TCells[celli] = mixture_.THE
        (
            hCells[celli],
            pCells[celli],
            TCells[celli]
        );

        rhoCells[celli] = volMixture_.rho(pCells[celli], TCells[celli]);

        alphaCells[celli] =
            volMixture_.kappa(pCells[celli], TCells[celli])
            /
            mixture_.Cpv(pCells[celli], TCells[celli]);
    }

    volScalarField::Boundary& pBf =
        this->p_.boundaryFieldRef();

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
        fvPatchScalarField& pp = pBf[patchi];
        fvPatchScalarField& pT = TBf[patchi];
        fvPatchScalarField& prho = rhoBf[patchi];
        fvPatchScalarField& phe = heBf[patchi];
        fvPatchScalarField& palpha = alphaBf[patchi];

        if (pT.fixesValue())
        {
            forAll(pT, facei)
            {
                const typename MixtureType::thermoType& mixture_ =
                    this->patchFaceMixture(patchi, facei);

                const typename MixtureType::thermoType& volMixture_ =
                    this->patchFaceVolMixture
                    (
                        pp[facei],
                        pT[facei],
                        patchi,
                        facei
                    );


                phe[facei] = mixture_.HE(pp[facei], pT[facei]);
                prho[facei] = volMixture_.rho(pp[facei], pT[facei]);

                palpha[facei] =
                    volMixture_.kappa(pp[facei], pT[facei])
                  / mixture_.Cpv(pp[facei], pT[facei]);
            }
        }
        else
        {
            forAll(pT, facei)
            {
                const typename MixtureType::thermoType& mixture_ =
                    this->patchFaceMixture(patchi, facei);

                const typename MixtureType::thermoType& volMixture_ =
                    this->patchFaceVolMixture
                    (
                        pp[facei],
                        pT[facei],
                        patchi,
                        facei
                    );

                pT[facei] = mixture_.THE(phe[facei], pp[facei] ,pT[facei]);
                prho[facei] = volMixture_.rho(pp[facei], pT[facei]);

                palpha[facei] =
                    volMixture_.kappa(pp[facei], pT[facei])
                  / mixture_.Cpv(pp[facei], pT[facei]);
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

    volVectorField& Kappa = tKappa.ref();
    vectorField& KappaCells = Kappa.primitiveFieldRef();
    const scalarField& TCells = this->T_;
    const scalarField& pCells = this->p_;

    forAll(KappaCells, celli)
    {
        Kappa[celli] =
            this->cellVolMixture
            (
                pCells[celli],
                TCells[celli],
                celli
            ).Kappa(pCells[celli], TCells[celli]);
    }

    volVectorField::Boundary& KappaBf = Kappa.boundaryFieldRef();

    forAll(KappaBf, patchi)
    {
        vectorField& Kappap = KappaBf[patchi];
        const scalarField& pT = this->T_.boundaryField()[patchi];
        const scalarField& pp = this->p_.boundaryField()[patchi];

        forAll(Kappap, facei)
        {
            Kappap[facei] =
                this->patchFaceVolMixture
                (
                    pp[facei],
                    pT[facei],
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
    const scalarField& pp = this->p_.boundaryField()[patchi];
    const scalarField& Tp = this->T_.boundaryField()[patchi];
    tmp<vectorField> tKappa(new vectorField(pp.size()));

    vectorField& Kappap = tKappa.ref();

    forAll(Tp, facei)
    {
        Kappap[facei] =
            this->patchFaceVolMixture
            (
                pp[facei],
                Tp[facei],
                patchi,
                facei
            ).Kappa(pp[facei], Tp[facei]);
    }

    return tKappa;
}


// ************************************************************************* //
