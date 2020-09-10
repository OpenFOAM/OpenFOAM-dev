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

#include "heThermo.H"
#include "gradientEnergyFvPatchScalarField.H"
#include "mixedEnergyFvPatchScalarField.H"

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class BasicThermo, class MixtureType>
template
<
    class CellMixture,
    class PatchFaceMixture,
    class Method,
    class ... Args
>
Foam::tmp<Foam::volScalarField>
Foam::heThermo<BasicThermo, MixtureType>::volScalarFieldProperty
(
    const word& psiName,
    const dimensionSet& psiDim,
    CellMixture cellMixture,
    PatchFaceMixture patchFaceMixture,
    Method psiMethod,
    const Args& ... args
) const
{
    tmp<volScalarField> tPsi
    (
        volScalarField::New
        (
            IOobject::groupName(psiName, this->group()),
            this->T_.mesh(),
            psiDim
        )
    );

    volScalarField& psi = tPsi.ref();

    forAll(this->T_, celli)
    {
        psi[celli] = ((this->*cellMixture)(celli).*psiMethod)(args[celli] ...);
    }

    volScalarField::Boundary& psiBf = psi.boundaryFieldRef();

    forAll(psiBf, patchi)
    {
        fvPatchScalarField& pPsi = psiBf[patchi];

        forAll(this->T_.boundaryField()[patchi], facei)
        {
            pPsi[facei] =
                ((this->*patchFaceMixture)(patchi, facei).*psiMethod)
                (
                    args.boundaryField()[patchi][facei] ...
                );
        }
    }

    return tPsi;
}


template<class BasicThermo, class MixtureType>
template<class CellMixture, class Method, class ... Args>
Foam::tmp<Foam::scalarField>
Foam::heThermo<BasicThermo, MixtureType>::cellSetProperty
(
    CellMixture cellMixture,
    Method psiMethod,
    const labelList& cells,
    const Args& ... args
) const
{
    // Note: Args are fields for the set, not for the mesh as a whole. The
    // cells list is only used to get the mixture.

    tmp<scalarField> tPsi(new scalarField(cells.size()));
    scalarField& psi = tPsi.ref();

    forAll(cells, celli)
    {
        psi[celli] =
            ((this->*cellMixture)(cells[celli]).*psiMethod)(args[celli] ...);
    }

    return tPsi;
}


template<class BasicThermo, class MixtureType>
template<class PatchFaceMixture, class Method, class ... Args>
Foam::tmp<Foam::scalarField>
Foam::heThermo<BasicThermo, MixtureType>::patchFieldProperty
(
    PatchFaceMixture patchFaceMixture,
    Method psiMethod,
    const label patchi,
    const Args& ... args
) const
{
    tmp<scalarField> tPsi
    (
        new scalarField(this->T_.boundaryField()[patchi].size())
    );
    scalarField& psi = tPsi.ref();

    forAll(this->T_.boundaryField()[patchi], facei)
    {
        psi[facei] =
            ((this->*patchFaceMixture)(patchi, facei).*psiMethod)
            (
                args[facei] ...
            );
    }

    return tPsi;
}


template<class BasicThermo, class MixtureType>
Foam::UIndirectList<Foam::scalar>
Foam::heThermo<BasicThermo, MixtureType>::cellSetScalarList
(
    const volScalarField& psi,
    const labelList& cells
)
{
    return UIndirectList<scalar>(psi, cells);
}


template<class BasicThermo, class MixtureType>
Foam::UniformField<Foam::scalar>
Foam::heThermo<BasicThermo, MixtureType>::cellSetScalarList
(
    const uniformGeometricScalarField& psi,
    const labelList&
)
{
    return psi.primitiveField();
}


template<class BasicThermo, class MixtureType>
void Foam::heThermo<BasicThermo, MixtureType>::
heBoundaryCorrection(volScalarField& h)
{
    volScalarField::Boundary& hBf = h.boundaryFieldRef();

    forAll(hBf, patchi)
    {
        if (isA<gradientEnergyFvPatchScalarField>(hBf[patchi]))
        {
            refCast<gradientEnergyFvPatchScalarField>(hBf[patchi]).gradient()
                = hBf[patchi].fvPatchField::snGrad();
        }
        else if (isA<mixedEnergyFvPatchScalarField>(hBf[patchi]))
        {
            refCast<mixedEnergyFvPatchScalarField>(hBf[patchi]).refGrad()
                = hBf[patchi].fvPatchField::snGrad();
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicThermo, class MixtureType>
Foam::heThermo<BasicThermo, MixtureType>::heThermo
(
    const fvMesh& mesh,
    const word& phaseName
)
:
    BasicThermo(mesh, phaseName),
    MixtureType(*this, mesh, phaseName),

    he_
    (
        IOobject
        (
            BasicThermo::phasePropertyName
            (
                MixtureType::thermoType::heName(),
                phaseName
            ),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        volScalarFieldProperty
        (
            "he",
            dimEnergy/dimMass,
            &MixtureType::cellThermoMixture,
            &MixtureType::patchFaceThermoMixture,
            &MixtureType::thermoMixtureType::HE,
            this->p_,
            this->T_
        ),
        this->heBoundaryTypes(),
        this->heBoundaryBaseTypes()
    )
{
    heBoundaryCorrection(he_);
}


template<class BasicThermo, class MixtureType>
Foam::heThermo<BasicThermo, MixtureType>::heThermo
(
    const fvMesh& mesh,
    const dictionary& dict,
    const word& phaseName
)
:
    BasicThermo(mesh, dict, phaseName),
    MixtureType(*this, mesh, phaseName),

    he_
    (
        IOobject
        (
            BasicThermo::phasePropertyName
            (
                MixtureType::thermoType::heName(),
                phaseName
            ),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        volScalarFieldProperty
        (
            "he",
            dimEnergy/dimMass,
            &MixtureType::cellThermoMixture,
            &MixtureType::patchFaceThermoMixture,
            &MixtureType::thermoMixtureType::HE,
            this->p_,
            this->T_
        ),
        this->heBoundaryTypes(),
        this->heBoundaryBaseTypes()
    )
{
    heBoundaryCorrection(he_);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class BasicThermo, class MixtureType>
Foam::heThermo<BasicThermo, MixtureType>::~heThermo()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicThermo, class MixtureType>
Foam::tmp<Foam::volScalarField> Foam::heThermo<BasicThermo, MixtureType>::he
(
    const volScalarField& p,
    const volScalarField& T
) const
{
    return volScalarFieldProperty
    (
        "he",
        dimEnergy/dimMass,
        &MixtureType::cellThermoMixture,
        &MixtureType::patchFaceThermoMixture,
        &MixtureType::thermoMixtureType::HE,
        p,
        T
    );
}


template<class BasicThermo, class MixtureType>
Foam::tmp<Foam::scalarField> Foam::heThermo<BasicThermo, MixtureType>::he
(
    const scalarField& T,
    const labelList& cells
) const
{
    return cellSetProperty
    (
        &MixtureType::cellThermoMixture,
        &MixtureType::thermoMixtureType::HE,
        cells,
        cellSetScalarList(this->p_, cells),
        T
    );
}


template<class BasicThermo, class MixtureType>
Foam::tmp<Foam::scalarField> Foam::heThermo<BasicThermo, MixtureType>::he
(
    const scalarField& T,
    const label patchi
) const
{
    return patchFieldProperty
    (
        &MixtureType::patchFaceThermoMixture,
        &MixtureType::thermoMixtureType::HE,
        patchi,
        this->p_.boundaryField()[patchi],
        T
    );
}


template<class BasicThermo, class MixtureType>
Foam::tmp<Foam::volScalarField>
Foam::heThermo<BasicThermo, MixtureType>::hs() const
{
    return volScalarFieldProperty
    (
        "hs",
        dimEnergy/dimMass,
        &MixtureType::cellThermoMixture,
        &MixtureType::patchFaceThermoMixture,
        &MixtureType::thermoMixtureType::Hs,
        this->p_,
        this->T_
    );
}


template<class BasicThermo, class MixtureType>
Foam::tmp<Foam::volScalarField> Foam::heThermo<BasicThermo, MixtureType>::hs
(
    const volScalarField& p,
    const volScalarField& T
) const
{
    return volScalarFieldProperty
    (
        "hs",
        dimEnergy/dimMass,
        &MixtureType::cellThermoMixture,
        &MixtureType::patchFaceThermoMixture,
        &MixtureType::thermoMixtureType::Hs,
        p,
        T
    );
}


template<class BasicThermo, class MixtureType>
Foam::tmp<Foam::scalarField> Foam::heThermo<BasicThermo, MixtureType>::hs
(
    const scalarField& T,
    const labelList& cells
) const
{
    return cellSetProperty
    (
        &MixtureType::cellThermoMixture,
        &MixtureType::thermoMixtureType::Hs,
        cells,
        cellSetScalarList(this->p_, cells),
        T
    );
}


template<class BasicThermo, class MixtureType>
Foam::tmp<Foam::scalarField> Foam::heThermo<BasicThermo, MixtureType>::hs
(
    const scalarField& T,
    const label patchi
) const
{
    return patchFieldProperty
    (
        &MixtureType::patchFaceThermoMixture,
        &MixtureType::thermoMixtureType::Hs,
        patchi,
        this->p_.boundaryField()[patchi],
        T
    );
}


template<class BasicThermo, class MixtureType>
Foam::tmp<Foam::volScalarField>
Foam::heThermo<BasicThermo, MixtureType>::ha() const
{
    return volScalarFieldProperty
    (
        "ha",
        dimEnergy/dimMass,
        &MixtureType::cellThermoMixture,
        &MixtureType::patchFaceThermoMixture,
        &MixtureType::thermoMixtureType::Ha,
        this->p_,
        this->T_
    );
}


template<class BasicThermo, class MixtureType>
Foam::tmp<Foam::volScalarField> Foam::heThermo<BasicThermo, MixtureType>::ha
(
    const volScalarField& p,
    const volScalarField& T
) const
{
    return volScalarFieldProperty
    (
        "ha",
        dimEnergy/dimMass,
        &MixtureType::cellThermoMixture,
        &MixtureType::patchFaceThermoMixture,
        &MixtureType::thermoMixtureType::Ha,
        p,
        T
    );
}


template<class BasicThermo, class MixtureType>
Foam::tmp<Foam::scalarField> Foam::heThermo<BasicThermo, MixtureType>::ha
(
    const scalarField& T,
    const labelList& cells
) const
{
    return cellSetProperty
    (
        &MixtureType::cellThermoMixture,
        &MixtureType::thermoMixtureType::Ha,
        cells,
        cellSetScalarList(this->p_, cells),
        T
    );
}


template<class BasicThermo, class MixtureType>
Foam::tmp<Foam::scalarField> Foam::heThermo<BasicThermo, MixtureType>::ha
(
    const scalarField& T,
    const label patchi
) const
{
    return patchFieldProperty
    (
        &MixtureType::patchFaceThermoMixture,
        &MixtureType::thermoMixtureType::Ha,
        patchi,
        this->p_.boundaryField()[patchi],
        T
    );
}


template<class BasicThermo, class MixtureType>
Foam::tmp<Foam::volScalarField>
Foam::heThermo<BasicThermo, MixtureType>::hc() const
{
    return volScalarFieldProperty
    (
        "hc",
        dimEnergy/dimMass,
        &MixtureType::cellThermoMixture,
        &MixtureType::patchFaceThermoMixture,
        &MixtureType::thermoMixtureType::Hf
    );
}


template<class BasicThermo, class MixtureType>
Foam::tmp<Foam::scalarField> Foam::heThermo<BasicThermo, MixtureType>::Cp
(
    const scalarField& T,
    const label patchi
) const
{
    return patchFieldProperty
    (
        &MixtureType::patchFaceThermoMixture,
        &MixtureType::thermoMixtureType::Cp,
        patchi,
        this->p_.boundaryField()[patchi],
        T
    );
}


template<class BasicThermo, class MixtureType>
Foam::tmp<Foam::volScalarField>
Foam::heThermo<BasicThermo, MixtureType>::Cp() const
{
    return volScalarFieldProperty
    (
        "Cp",
        dimEnergy/dimMass/dimTemperature,
        &MixtureType::cellThermoMixture,
        &MixtureType::patchFaceThermoMixture,
        &MixtureType::thermoMixtureType::Cp,
        this->p_,
        this->T_
    );
}


template<class BasicThermo, class MixtureType>
Foam::tmp<Foam::scalarField>
Foam::heThermo<BasicThermo, MixtureType>::Cv
(
    const scalarField& T,
    const label patchi
) const
{
    return patchFieldProperty
    (
        &MixtureType::patchFaceThermoMixture,
        &MixtureType::thermoMixtureType::Cv,
        patchi,
        this->p_.boundaryField()[patchi],
        T
    );
}


template<class BasicThermo, class MixtureType>
Foam::tmp<Foam::volScalarField>
Foam::heThermo<BasicThermo, MixtureType>::Cv() const
{
    return volScalarFieldProperty
    (
        "Cv",
        dimEnergy/dimMass/dimTemperature,
        &MixtureType::cellThermoMixture,
        &MixtureType::patchFaceThermoMixture,
        &MixtureType::thermoMixtureType::Cv,
        this->p_,
        this->T_
    );
}


template<class BasicThermo, class MixtureType>
Foam::tmp<Foam::scalarField> Foam::heThermo<BasicThermo, MixtureType>::gamma
(
    const scalarField& T,
    const label patchi
) const
{
    return patchFieldProperty
    (
        &MixtureType::patchFaceThermoMixture,
        &MixtureType::thermoMixtureType::gamma,
        patchi,
        this->p_.boundaryField()[patchi],
        T
    );
}


template<class BasicThermo, class MixtureType>
Foam::tmp<Foam::volScalarField>
Foam::heThermo<BasicThermo, MixtureType>::gamma() const
{
    return volScalarFieldProperty
    (
        "gamma",
        dimless,
        &MixtureType::cellThermoMixture,
        &MixtureType::patchFaceThermoMixture,
        &MixtureType::thermoMixtureType::gamma,
        this->p_,
        this->T_
    );
}


template<class BasicThermo, class MixtureType>
Foam::tmp<Foam::scalarField> Foam::heThermo<BasicThermo, MixtureType>::Cpv
(
    const scalarField& T,
    const label patchi
) const
{
    return patchFieldProperty
    (
        &MixtureType::patchFaceThermoMixture,
        &MixtureType::thermoMixtureType::Cpv,
        patchi,
        this->p_.boundaryField()[patchi],
        T
    );
}


template<class BasicThermo, class MixtureType>
Foam::tmp<Foam::volScalarField>
Foam::heThermo<BasicThermo, MixtureType>::Cpv() const
{
    return volScalarFieldProperty
    (
        "Cpv",
        dimEnergy/dimMass/dimTemperature,
        &MixtureType::cellThermoMixture,
        &MixtureType::patchFaceThermoMixture,
        &MixtureType::thermoMixtureType::Cpv,
        this->p_,
        this->T_
    );
}


template<class BasicThermo, class MixtureType>
Foam::tmp<Foam::scalarField> Foam::heThermo<BasicThermo, MixtureType>::CpByCpv
(
    const scalarField& T,
    const label patchi
) const
{
    return patchFieldProperty
    (
        &MixtureType::patchFaceThermoMixture,
        &MixtureType::thermoMixtureType::CpByCpv,
        patchi,
        this->p_.boundaryField()[patchi],
        T
    );
}


template<class BasicThermo, class MixtureType>
Foam::tmp<Foam::volScalarField>
Foam::heThermo<BasicThermo, MixtureType>::CpByCpv() const
{
    return volScalarFieldProperty
    (
        "CpByCpv",
        dimless,
        &MixtureType::cellThermoMixture,
        &MixtureType::patchFaceThermoMixture,
        &MixtureType::thermoMixtureType::CpByCpv,
        this->p_,
        this->T_
    );
}


template<class BasicThermo, class MixtureType>
Foam::tmp<Foam::volScalarField> Foam::heThermo<BasicThermo, MixtureType>::THE
(
    const volScalarField& h,
    const volScalarField& p,
    const volScalarField& T0
) const
{
    return volScalarFieldProperty
    (
        "T",
        dimTemperature,
        &MixtureType::cellThermoMixture,
        &MixtureType::patchFaceThermoMixture,
        &MixtureType::thermoMixtureType::THE,
        h,
        p,
        T0
    );
}


template<class BasicThermo, class MixtureType>
Foam::tmp<Foam::scalarField> Foam::heThermo<BasicThermo, MixtureType>::THE
(
    const scalarField& h,
    const scalarField& T0,
    const labelList& cells
) const
{
    return cellSetProperty
    (
        &MixtureType::cellThermoMixture,
        &MixtureType::thermoMixtureType::THE,
        cells,
        h,
        cellSetScalarList(this->p_, cells),
        T0
    );
}


template<class BasicThermo, class MixtureType>
Foam::tmp<Foam::scalarField> Foam::heThermo<BasicThermo, MixtureType>::THE
(
    const scalarField& h,
    const scalarField& T0,
    const label patchi
) const
{
    return patchFieldProperty
    (
        &MixtureType::patchFaceThermoMixture,
        &MixtureType::thermoMixtureType::THE,
        patchi,
        h,
        this->p_.boundaryField()[patchi],
        T0
    );
}


template<class BasicThermo, class MixtureType>
Foam::tmp<Foam::volScalarField>
Foam::heThermo<BasicThermo, MixtureType>::W() const
{
    return volScalarFieldProperty
    (
        "W",
        dimMass/dimMoles,
        &MixtureType::cellThermoMixture,
        &MixtureType::patchFaceThermoMixture,
        &MixtureType::thermoMixtureType::W
    );
}


template<class BasicThermo, class MixtureType>
Foam::tmp<Foam::scalarField> Foam::heThermo<BasicThermo, MixtureType>::W
(
    const label patchi
) const
{
    return patchFieldProperty
    (
        &MixtureType::patchFaceThermoMixture,
        &MixtureType::thermoMixtureType::W,
        patchi
    );
}


template<class BasicThermo, class MixtureType>
Foam::tmp<Foam::volScalarField>
Foam::heThermo<BasicThermo, MixtureType>::kappa() const
{
    return volScalarField::New
    (
        "kappa",
        Cp()*this->alpha_
    );
}


template<class BasicThermo, class MixtureType>
Foam::tmp<Foam::scalarField> Foam::heThermo<BasicThermo, MixtureType>::kappa
(
    const label patchi
) const
{
    return
        Cp
        (
            this->T_.boundaryField()[patchi],
            patchi
        )*this->alpha_.boundaryField()[patchi];
}


template<class BasicThermo, class MixtureType>
Foam::tmp<Foam::volScalarField>
Foam::heThermo<BasicThermo, MixtureType>::alphahe() const
{
    return volScalarField::New
    (
        "alphahe",
        this->CpByCpv()*this->alpha_
    );
}


template<class BasicThermo, class MixtureType>
Foam::tmp<Foam::scalarField>
Foam::heThermo<BasicThermo, MixtureType>::alphahe(const label patchi) const
{
    return
        this->CpByCpv
        (
            this->T_.boundaryField()[patchi],
            patchi
        )
       *this->alpha_.boundaryField()[patchi];
}


template<class BasicThermo, class MixtureType>
Foam::tmp<Foam::volScalarField>
Foam::heThermo<BasicThermo, MixtureType>::kappaEff
(
    const volScalarField& alphat
) const
{
    return volScalarField::New
    (
        "kappaEff",
        Cp()*(this->alpha_ + alphat)
    );
}


template<class BasicThermo, class MixtureType>
Foam::tmp<Foam::scalarField>
Foam::heThermo<BasicThermo, MixtureType>::kappaEff
(
    const scalarField& alphat,
    const label patchi
) const
{
    return
        Cp
        (
            this->T_.boundaryField()[patchi],
            patchi
        )
       *(this->alpha_.boundaryField()[patchi] + alphat);
}


template<class BasicThermo, class MixtureType>
Foam::tmp<Foam::volScalarField>
Foam::heThermo<BasicThermo, MixtureType>::alphaEff
(
    const volScalarField& alphat
) const
{
    return volScalarField::New
    (
        "alphaEff",
        this->CpByCpv()*(this->alpha_ + alphat)
    );
}


template<class BasicThermo, class MixtureType>
Foam::tmp<Foam::scalarField>
Foam::heThermo<BasicThermo, MixtureType>::alphaEff
(
    const scalarField& alphat,
    const label patchi
) const
{
    return
        this->CpByCpv
        (
            this->T_.boundaryField()[patchi],
            patchi
        )
       *(this->alpha_.boundaryField()[patchi] + alphat);
}


template<class BasicThermo, class MixtureType>
bool Foam::heThermo<BasicThermo, MixtureType>::read()
{
    if (BasicThermo::read())
    {
        MixtureType::read(*this);
        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
