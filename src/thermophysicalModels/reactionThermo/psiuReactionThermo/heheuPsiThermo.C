/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
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

#include "heheuPsiThermo.H"
#include "fvMesh.H"
#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class BasicPsiThermo, class MixtureType>
void Foam::heheuPsiThermo<BasicPsiThermo, MixtureType>::calculate()
{
    const scalarField& hCells = this->he_;
    const scalarField& heuCells = this->heu_;
    const scalarField& pCells = this->p_;

    scalarField& TCells = this->T_.primitiveFieldRef();
    scalarField& TuCells = this->Tu_.primitiveFieldRef();
    scalarField& psiCells = this->psi_.primitiveFieldRef();
    scalarField& muCells = this->mu_.primitiveFieldRef();
    scalarField& alphaCells = this->alpha_.primitiveFieldRef();

    forAll(TCells, celli)
    {
        const typename MixtureType::thermoType& mixture_ =
            this->cellMixture(celli);

        TCells[celli] = mixture_.THE
        (
            hCells[celli],
            pCells[celli],
            TCells[celli]
        );

        psiCells[celli] = mixture_.psi(pCells[celli], TCells[celli]);

        muCells[celli] = mixture_.mu(pCells[celli], TCells[celli]);
        alphaCells[celli] = mixture_.alphah(pCells[celli], TCells[celli]);

        TuCells[celli] = this->cellReactants(celli).THE
        (
            heuCells[celli],
            pCells[celli],
            TuCells[celli]
        );
    }

    volScalarField::Boundary& pBf =
        this->p_.boundaryFieldRef();

    volScalarField::Boundary& TBf =
        this->T_.boundaryFieldRef();

    volScalarField::Boundary& TuBf =
        this->Tu_.boundaryFieldRef();

    volScalarField::Boundary& psiBf =
        this->psi_.boundaryFieldRef();

    volScalarField::Boundary& heBf =
        this->he().boundaryFieldRef();

    volScalarField::Boundary& heuBf =
        this->heu().boundaryFieldRef();

    volScalarField::Boundary& muBf =
        this->mu_.boundaryFieldRef();

    volScalarField::Boundary& alphaBf =
        this->alpha_.boundaryFieldRef();

    forAll(this->T_.boundaryField(), patchi)
    {
        fvPatchScalarField& pp = pBf[patchi];
        fvPatchScalarField& pT = TBf[patchi];
        fvPatchScalarField& pTu = TuBf[patchi];
        fvPatchScalarField& ppsi = psiBf[patchi];
        fvPatchScalarField& phe = heBf[patchi];
        fvPatchScalarField& pheu = heuBf[patchi];
        fvPatchScalarField& pmu = muBf[patchi];
        fvPatchScalarField& palpha = alphaBf[patchi];

        if (pT.fixesValue())
        {
            forAll(pT, facei)
            {
                const typename MixtureType::thermoType& mixture_ =
                    this->patchFaceMixture(patchi, facei);

                phe[facei] = mixture_.HE(pp[facei], pT[facei]);

                ppsi[facei] = mixture_.psi(pp[facei], pT[facei]);
                pmu[facei] = mixture_.mu(pp[facei], pT[facei]);
                palpha[facei] = mixture_.alphah(pp[facei], pT[facei]);
            }
        }
        else
        {
            forAll(pT, facei)
            {
                const typename MixtureType::thermoType& mixture_ =
                    this->patchFaceMixture(patchi, facei);

                pT[facei] = mixture_.THE(phe[facei], pp[facei], pT[facei]);

                ppsi[facei] = mixture_.psi(pp[facei], pT[facei]);
                pmu[facei] = mixture_.mu(pp[facei], pT[facei]);
                palpha[facei] = mixture_.alphah(pp[facei], pT[facei]);

                pTu[facei] =
                    this->patchFaceReactants(patchi, facei)
                   .THE(pheu[facei], pp[facei], pTu[facei]);
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicPsiThermo, class MixtureType>
Foam::heheuPsiThermo<BasicPsiThermo, MixtureType>::heheuPsiThermo
(
    const fvMesh& mesh,
    const word& phaseName
)
:
    heThermo<psiuReactionThermo, MixtureType>(mesh, phaseName),
    Tu_
    (
        IOobject
        (
            "Tu",
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ),
    heu_
    (
        IOobject
        (
            MixtureType::thermoType::heName() + 'u',
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        this->volScalarFieldProperty
        (
            MixtureType::thermoType::heName() + 'u',
            dimEnergy/dimMass,
            &MixtureType::cellReactants,
            &MixtureType::patchFaceReactants,
            &MixtureType::thermoType::HE,
            this->p_,
            this->Tu_
        ),
        this->heuBoundaryTypes()
    )
{
    this->heuBoundaryCorrection(this->heu_);

    calculate();

    this->psi_.oldTime(); // Switch on saving old time
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class BasicPsiThermo, class MixtureType>
Foam::heheuPsiThermo<BasicPsiThermo, MixtureType>::~heheuPsiThermo()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicPsiThermo, class MixtureType>
void Foam::heheuPsiThermo<BasicPsiThermo, MixtureType>::correct()
{
    if (debug)
    {
        InfoInFunction << endl;
    }

    // force the saving of the old-time values
    this->psi_.oldTime();

    calculate();

    if (debug)
    {
        Info<< "    Finished" << endl;
    }
}


template<class BasicPsiThermo, class MixtureType>
Foam::tmp<Foam::scalarField>
Foam::heheuPsiThermo<BasicPsiThermo, MixtureType>::heu
(
    const scalarField& Tu,
    const labelList& cells
) const
{
    return this->cellSetProperty
    (
        &MixtureType::cellReactants,
        &MixtureType::thermoType::HE,
        cells,
        UIndirectList<scalar>(this->p(), cells),
        Tu
    );
}


template<class BasicPsiThermo, class MixtureType>
Foam::tmp<Foam::scalarField>
Foam::heheuPsiThermo<BasicPsiThermo, MixtureType>::heu
(
    const scalarField& Tu,
    const label patchi
) const
{
    return this->patchFieldProperty
    (
        &MixtureType::patchFaceReactants,
        &MixtureType::thermoType::HE,
        patchi,
        this->p().boundaryField()[patchi],
        Tu
    );
}


template<class BasicPsiThermo, class MixtureType>
Foam::tmp<Foam::volScalarField>
Foam::heheuPsiThermo<BasicPsiThermo, MixtureType>::Tb() const
{
    return this->volScalarFieldProperty
    (
        "Tb",
        dimTemperature,
        &MixtureType::cellProducts,
        &MixtureType::patchFaceProducts,
        &MixtureType::thermoType::THE,
        this->he_,
        this->p_,
        this->T_
    );
}


template<class BasicPsiThermo, class MixtureType>
Foam::tmp<Foam::volScalarField>
Foam::heheuPsiThermo<BasicPsiThermo, MixtureType>::psiu() const
{
    return this->volScalarFieldProperty
    (
        "psiu",
        this->psi_.dimensions(),
        &MixtureType::cellReactants,
        &MixtureType::patchFaceReactants,
        &MixtureType::thermoType::psi,
        this->p_,
        this->T_
    );
}


template<class BasicPsiThermo, class MixtureType>
Foam::tmp<Foam::volScalarField>
Foam::heheuPsiThermo<BasicPsiThermo, MixtureType>::psib() const
{
    return this->volScalarFieldProperty
    (
        "psib",
        this->psi_.dimensions(),
        &MixtureType::cellProducts,
        &MixtureType::patchFaceProducts,
        &MixtureType::thermoType::psi,
        this->p_,
        this->T_
    );
}


template<class BasicPsiThermo, class MixtureType>
Foam::tmp<Foam::volScalarField>
Foam::heheuPsiThermo<BasicPsiThermo, MixtureType>::muu() const
{
    return this->volScalarFieldProperty
    (
        "muu",
        dimDynamicViscosity,
        &MixtureType::cellReactants,
        &MixtureType::patchFaceReactants,
        &MixtureType::thermoType::mu,
        this->p_,
        this->T_
    );
}


template<class BasicPsiThermo, class MixtureType>
Foam::tmp<Foam::volScalarField>
Foam::heheuPsiThermo<BasicPsiThermo, MixtureType>::mub() const
{
    return this->volScalarFieldProperty
    (
        "mub",
        dimDynamicViscosity,
        &MixtureType::cellProducts,
        &MixtureType::patchFaceProducts,
        &MixtureType::thermoType::mu,
        this->p_,
        this->T_
    );
}


// ************************************************************************* //
