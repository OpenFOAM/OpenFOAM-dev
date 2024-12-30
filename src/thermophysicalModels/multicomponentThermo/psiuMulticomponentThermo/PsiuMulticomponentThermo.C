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

#include "PsiuMulticomponentThermo.H"
#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class BaseThermo>
void Foam::PsiuMulticomponentThermo<BaseThermo>::calculate()
{
    const scalarField& hCells = this->he_;
    const scalarField& heuCells = this->heu_;
    const scalarField& pCells = this->p_;

    scalarField& TCells = this->T_.primitiveFieldRef();
    scalarField& TuCells = this->Tu_.primitiveFieldRef();
    scalarField& CpCells = this->Cp_.primitiveFieldRef();
    scalarField& CvCells = this->Cv_.primitiveFieldRef();
    scalarField& psiCells = this->psi_.primitiveFieldRef();
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

        muCells[celli] = transportMixture.mu(pCells[celli], TCells[celli]);
        kappaCells[celli] =
            transportMixture.kappa(pCells[celli], TCells[celli]);

        TuCells[celli] = this->reactants(composition).The
        (
            heuCells[celli],
            pCells[celli],
            TuCells[celli]
        );
    }

    volScalarField::Boundary& pBf = this->p_.boundaryFieldRef();
    volScalarField::Boundary& TBf = this->T_.boundaryFieldRef();
    volScalarField::Boundary& TuBf = this->Tu_.boundaryFieldRef();
    volScalarField::Boundary& CpBf = this->Cp_.boundaryFieldRef();
    volScalarField::Boundary& CvBf = this->Cv_.boundaryFieldRef();
    volScalarField::Boundary& psiBf = this->psi_.boundaryFieldRef();
    volScalarField::Boundary& heBf = this->he().boundaryFieldRef();
    volScalarField::Boundary& heuBf = this->heu().boundaryFieldRef();
    volScalarField::Boundary& muBf = this->mu_.boundaryFieldRef();
    volScalarField::Boundary& kappaBf = this->kappa_.boundaryFieldRef();

    forAll(TBf, patchi)
    {
        fvPatchScalarField& pPf = pBf[patchi];
        fvPatchScalarField& TPf = TBf[patchi];
        fvPatchScalarField& TuPf = TuBf[patchi];
        fvPatchScalarField& CpPf = CpBf[patchi];
        fvPatchScalarField& CvPf = CvBf[patchi];
        fvPatchScalarField& psiPf = psiBf[patchi];
        fvPatchScalarField& hePf = heBf[patchi];
        fvPatchScalarField& heuPf = heuBf[patchi];
        fvPatchScalarField& muPf = muBf[patchi];
        fvPatchScalarField& kappaPf = kappaBf[patchi];

        if (TPf.fixesValue())
        {
            forAll(TPf, facei)
            {
                auto composition =
                    this->patchFaceComposition(Yslicer, patchi, facei);

                const typename BaseThermo::mixtureType::thermoMixtureType&
                    thermoMixture = this->thermoMixture(composition);

                const typename BaseThermo::mixtureType::transportMixtureType&
                    transportMixture =
                    this->transportMixture(composition, thermoMixture);

                hePf[facei] = thermoMixture.he(pPf[facei], TPf[facei]);

                CpPf[facei] = thermoMixture.Cp(pPf[facei], TPf[facei]);
                CvPf[facei] = thermoMixture.Cv(pPf[facei], TPf[facei]);
                psiPf[facei] = thermoMixture.psi(pPf[facei], TPf[facei]);
                muPf[facei] = transportMixture.mu(pPf[facei], TPf[facei]);
                kappaPf[facei] = transportMixture.kappa(pPf[facei], TPf[facei]);
            }
        }
        else
        {
            forAll(TPf, facei)
            {
                auto composition =
                    this->patchFaceComposition(Yslicer, patchi, facei);

                const typename BaseThermo::mixtureType::thermoMixtureType&
                    thermoMixture = this->thermoMixture(composition);

                const typename BaseThermo::mixtureType::transportMixtureType&
                    transportMixture =
                    this->transportMixture(composition, thermoMixture);

                TPf[facei] =
                    thermoMixture.The(hePf[facei], pPf[facei], TPf[facei]);

                CpPf[facei] = thermoMixture.Cp(pPf[facei], TPf[facei]);
                CvPf[facei] = thermoMixture.Cv(pPf[facei], TPf[facei]);
                psiPf[facei] = thermoMixture.psi(pPf[facei], TPf[facei]);
                muPf[facei] = transportMixture.mu(pPf[facei], TPf[facei]);
                kappaPf[facei] = transportMixture.kappa(pPf[facei], TPf[facei]);

                TuPf[facei] =
                    this->reactants(composition)
                   .The(heuPf[facei], pPf[facei], TuPf[facei]);
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BaseThermo>
Foam::PsiuMulticomponentThermo<BaseThermo>::PsiuMulticomponentThermo
(
    const fvMesh& mesh,
    const word& phaseName
)
:
    BaseThermo(mesh, phaseName),
    Tu_
    (
        IOobject
        (
            "Tu",
            mesh.time().name(),
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
            BaseThermo::mixtureType::thermoType::heName() + 'u',
            mesh.time().name(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        this->volScalarFieldProperty
        (
            BaseThermo::mixtureType::thermoType::heName() + 'u',
            dimEnergy/dimMass,
            &BaseThermo::mixtureType::reactants,
            &BaseThermo::mixtureType::thermoMixtureType::he,
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

template<class BaseThermo>
Foam::PsiuMulticomponentThermo<BaseThermo>::~PsiuMulticomponentThermo()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BaseThermo>
void Foam::PsiuMulticomponentThermo<BaseThermo>::correct()
{
    if (BaseThermo::debug)
    {
        InfoInFunction << endl;
    }

    // force the saving of the old-time values
    this->psi_.oldTime();

    calculate();

    if (BaseThermo::debug)
    {
        Info<< "    Finished" << endl;
    }
}


template<class BaseThermo>
Foam::tmp<Foam::volScalarField>
Foam::PsiuMulticomponentThermo<BaseThermo>::fres() const
{
    tmp<volScalarField> tfres
    (
        volScalarField::New
        (
            "fres",
            this->mesh(),
            dimless
        )
    );

    auto Yslicer = this->Yslicer();

    scalarField& fresCells = tfres.ref().primitiveFieldRef();

    forAll(fresCells, celli)
    {
        fresCells[celli] = BaseThermo::mixtureType::fres
        (
            this->cellComposition(Yslicer, celli)
        );
    }

    volScalarField::Boundary& fresBf = tfres.ref().boundaryFieldRef();

    forAll(fresBf, patchi)
    {
        fvPatchScalarField& fresPf = fresBf[patchi];

        forAll(fresPf, facei)
        {
            fresPf[facei] = BaseThermo::mixtureType::fres
            (
                this->patchFaceComposition(Yslicer, patchi, facei)
            );
        }
    }

    return tfres;
}


template<class BaseThermo>
void Foam::PsiuMulticomponentThermo<BaseThermo>::reset()
{
    BaseThermo::mixtureType::reset(this->Y());
}


template<class BaseThermo>
Foam::tmp<Foam::scalarField>
Foam::PsiuMulticomponentThermo<BaseThermo>::heu
(
    const scalarField& Tu,
    const labelList& cells
) const
{
    return this->cellSetProperty
    (
        &BaseThermo::mixtureType::reactants,
        &BaseThermo::mixtureType::thermoMixtureType::he,
        cells,
        UIndirectList<scalar>(this->p_, cells),
        Tu
    );
}


template<class BaseThermo>
Foam::tmp<Foam::scalarField>
Foam::PsiuMulticomponentThermo<BaseThermo>::heu
(
    const scalarField& Tu,
    const label patchi
) const
{
    return this->patchFieldProperty
    (
        &BaseThermo::mixtureType::reactants,
        &BaseThermo::mixtureType::thermoMixtureType::he,
        patchi,
        this->p_.boundaryField()[patchi],
        Tu
    );
}


template<class BaseThermo>
Foam::tmp<Foam::volScalarField>
Foam::PsiuMulticomponentThermo<BaseThermo>::Tb() const
{
    return this->volScalarFieldProperty
    (
        "Tb",
        dimTemperature,
        &BaseThermo::mixtureType::products,
        &BaseThermo::mixtureType::thermoMixtureType::The,
        this->he_,
        this->p_,
        this->T_
    );
}


template<class BaseThermo>
Foam::tmp<Foam::volScalarField>
Foam::PsiuMulticomponentThermo<BaseThermo>::hr() const
{
    return
        this->volScalarFieldProperty
        (
            "hf",
            dimEnergy/dimMass,
            &BaseThermo::mixtureType::reactants,
            &BaseThermo::mixtureType::thermoMixtureType::hf
        )
      - this->volScalarFieldProperty
        (
            "hf",
            dimEnergy/dimMass,
            &BaseThermo::mixtureType::products,
            &BaseThermo::mixtureType::thermoMixtureType::hf
        );
}


template<class BaseThermo>
Foam::tmp<Foam::volScalarField>
Foam::PsiuMulticomponentThermo<BaseThermo>::psiu() const
{
    return this->volScalarFieldProperty
    (
        "psiu",
        this->psi_.dimensions(),
        &BaseThermo::mixtureType::reactants,
        &BaseThermo::mixtureType::thermoMixtureType::psi,
        this->p_,
        this->Tu_
    );
}


template<class BaseThermo>
Foam::tmp<Foam::volScalarField>
Foam::PsiuMulticomponentThermo<BaseThermo>::psib() const
{
    const volScalarField Tb(this->Tb());

    return this->volScalarFieldProperty
    (
        "psib",
        this->psi_.dimensions(),
        &BaseThermo::mixtureType::products,
        &BaseThermo::mixtureType::thermoMixtureType::psi,
        this->p_,
        Tb
    );
}


template<class BaseThermo>
Foam::tmp<Foam::volScalarField>
Foam::PsiuMulticomponentThermo<BaseThermo>::muu() const
{
    return this->volScalarFieldProperty
    (
        "muu",
        dimDynamicViscosity,
        &BaseThermo::mixtureType::reactants,
        &BaseThermo::mixtureType::transportMixtureType::mu,
        this->p_,
        this->Tu_
    );
}


template<class BaseThermo>
Foam::tmp<Foam::volScalarField>
Foam::PsiuMulticomponentThermo<BaseThermo>::mub() const
{
    const volScalarField Tb(this->Tb());

    return this->volScalarFieldProperty
    (
        "mub",
        dimDynamicViscosity,
        &BaseThermo::mixtureType::products,
        &BaseThermo::mixtureType::transportMixtureType::mu,
        this->p_,
        Tb
    );
}


// ************************************************************************* //
