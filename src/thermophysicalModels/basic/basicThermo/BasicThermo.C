/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2025 OpenFOAM Foundation
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

#include "BasicThermo.H"
#include "gradientEnergyFvPatchScalarField.H"
#include "mixedEnergyFvPatchScalarField.H"

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class MixtureType, class BasicThermoType>
template<class Mixture, class Method, class ... Args>
Foam::tmp<Foam::volScalarField>
Foam::BasicThermo<MixtureType, BasicThermoType>::volScalarFieldProperty
(
    const word& psiName,
    const dimensionSet& psiDim,
    Mixture mixture,
    Method psiMethod,
    const Args& ... args
) const
{
    tmp<volScalarField> tPsi
    (
        volScalarField::New
        (
            IOobject::groupName(psiName, this->group()),
            this->mesh(),
            psiDim
        )
    );

    volScalarField& psi = tPsi.ref();

    auto Yslicer = this->Yslicer();

    forAll(psi, celli)
    {
        auto composition = this->cellComposition(Yslicer, celli);

        psi[celli] =
            ((this->*mixture)(composition).*psiMethod)(args[celli] ...);
    }

    volScalarField::Boundary& psiBf = psi.boundaryFieldRef();

    forAll(psiBf, patchi)
    {
        forAll(psiBf[patchi], patchFacei)
        {
            auto composition =
                this->patchFaceComposition(Yslicer, patchi, patchFacei);

            psiBf[patchi][patchFacei] =
                ((this->*mixture)(composition).*psiMethod)
                (
                    args.boundaryField()[patchi][patchFacei] ...
                );
        }
    }

    return tPsi;
}


template<class MixtureType, class BasicThermoType>
template<class Mixture, class Method, class ... Args>
Foam::tmp<Foam::volScalarField::Internal>
Foam::BasicThermo<MixtureType, BasicThermoType>::volInternalScalarFieldProperty
(
    const word& psiName,
    const dimensionSet& psiDim,
    Mixture mixture,
    Method psiMethod,
    const Args& ... args
) const
{
    tmp<volScalarField::Internal> tPsi
    (
        volScalarField::Internal::New
        (
            IOobject::groupName(psiName, this->group()),
            this->mesh(),
            psiDim
        )
    );

    volScalarField::Internal& psi = tPsi.ref();

    auto Yslicer = this->Yslicer();

    forAll(psi, celli)
    {
        auto composition = this->cellComposition(Yslicer, celli);

        psi[celli] =
            ((this->*mixture)(composition).*psiMethod)(args[celli] ...);
    }

    return tPsi;
}


template<class MixtureType, class BasicThermoType>
template<class Mixture, class Method, class ... Args>
Foam::tmp<Foam::scalarField>
Foam::BasicThermo<MixtureType, BasicThermoType>::cellSetProperty
(
    Mixture mixture,
    Method psiMethod,
    const labelList& cells,
    const Args& ... args
) const
{
    // Note: Args are fields for the set, not for the mesh as a whole. The
    // cells list is only used to get the mixture.

    tmp<scalarField> tPsi(new scalarField(cells.size()));
    scalarField& psi = tPsi.ref();

    auto Yslicer = this->Yslicer();

    forAll(cells, i)
    {
        auto composition = this->cellComposition(Yslicer, cells[i]);

        psi[i] = ((this->*mixture)(composition).*psiMethod)(args[i] ...);
    }

    return tPsi;
}


template<class MixtureType, class BasicThermoType>
template<class Mixture, class Method, class ... Args>
Foam::tmp<Foam::scalarField>
Foam::BasicThermo<MixtureType, BasicThermoType>::patchFieldProperty
(
    Mixture mixture,
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

    auto Yslicer = this->Yslicer();

    forAll(psi, patchFacei)
    {
        auto composition =
            this->patchFaceComposition(Yslicer, patchi, patchFacei);

        psi[patchFacei] =
            ((this->*mixture)(composition).*psiMethod)(args[patchFacei] ...);
    }

    return tPsi;
}


template<class MixtureType, class BasicThermoType>
template<class Mixture, class Method, class ... Args>
Foam::tmp<Foam::scalarField>
Foam::BasicThermo<MixtureType, BasicThermoType>::fieldSourceProperty
(
    Mixture mixture,
    Method psiMethod,
    const fvSource& source,
    const Args& ... args
) const
{
    const labelUList cells = source.cells();

    tmp<scalarField> tPsi(new scalarField(cells.size()));
    scalarField& psi = tPsi.ref();

    auto Yslicer = this->Yslicer(source);

    forAll(cells, i)
    {
        auto composition =
            this->sourceCellComposition(Yslicer, i);

        psi[i] =
            ((this->*mixture)(composition).*psiMethod)(args[i] ...);
    }

    return tPsi;
}


template<class MixtureType, class BasicThermoType>
Foam::UIndirectList<Foam::scalar>
Foam::BasicThermo<MixtureType, BasicThermoType>::cellSetScalarList
(
    const volScalarField& psi,
    const labelUList& cells
)
{
    return UIndirectList<scalar>(psi, cells);
}


template<class MixtureType, class BasicThermoType>
Foam::UniformField<Foam::scalar>
Foam::BasicThermo<MixtureType, BasicThermoType>::cellSetScalarList
(
    const uniformGeometricScalarField& psi,
    const labelUList&
)
{
    return psi.primitiveField();
}


template<class MixtureType, class BasicThermoType>
void Foam::BasicThermo<MixtureType, BasicThermoType>::
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

template<class MixtureType, class BasicThermoType>
Foam::BasicThermo<MixtureType, BasicThermoType>::BasicThermo
(
    const fvMesh& mesh,
    const word& phaseName
)
:
    physicalProperties(mesh, phaseName),
    MixtureType(properties()),
    BasicThermoType(properties(), mixture(), mesh, phaseName),

    he_
    (
        IOobject
        (
            BasicThermoType::phasePropertyName
            (
                MixtureType::thermoType::heName(),
                phaseName
            ),
            mesh.time().name(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        volScalarFieldProperty
        (
            "he",
            dimEnergy/dimMass,
            &MixtureType::thermoMixture,
            &MixtureType::thermoMixtureType::he,
            this->p_,
            this->T_
        ),
        this->heBoundaryTypes(),
        this->heBoundaryBaseTypes(),
        this->heSourcesTypes()
    ),

    Cp_
    (
        IOobject
        (
            BasicThermoType::phasePropertyName("Cp", phaseName),
            mesh.time().name(),
            mesh
        ),
        mesh,
        dimensionedScalar(dimEnergy/dimMass/dimTemperature, Zero)
    ),

    Cv_
    (
        IOobject
        (
            BasicThermoType::phasePropertyName("Cv", phaseName),
            mesh.time().name(),
            mesh
        ),
        mesh,
        dimensionedScalar(dimEnergy/dimMass/dimTemperature, Zero)
    )
{
    heBoundaryCorrection(he_);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class MixtureType, class BasicThermoType>
Foam::BasicThermo<MixtureType, BasicThermoType>::~BasicThermo()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class MixtureType, class BasicThermoType>
Foam::tmp<Foam::volScalarField>
Foam::BasicThermo<MixtureType, BasicThermoType>::W() const
{
    return volScalarFieldProperty
    (
        "W",
        dimMass/dimMoles,
        &MixtureType::thermoMixture,
        &MixtureType::thermoMixtureType::W
    );
}


template<class MixtureType, class BasicThermoType>
Foam::tmp<Foam::scalarField>
Foam::BasicThermo<MixtureType, BasicThermoType>::W
(
    const label patchi
) const
{
    return patchFieldProperty
    (
        &MixtureType::thermoMixture,
        &MixtureType::thermoMixtureType::W,
        patchi
    );
}


template<class MixtureType, class BasicThermoType>
const Foam::volScalarField&
Foam::BasicThermo<MixtureType, BasicThermoType>::Cpv() const
{
    if (MixtureType::thermoType::enthalpy())
    {
        return Cp_;
    }
    else
    {
        return Cv_;
    }
}


template<class MixtureType, class BasicThermoType>
Foam::tmp<Foam::volScalarField>
Foam::BasicThermo<MixtureType, BasicThermoType>::he
(
    const volScalarField& p,
    const volScalarField& T
) const
{
    return volScalarFieldProperty
    (
        "he",
        dimEnergy/dimMass,
        &MixtureType::thermoMixture,
        &MixtureType::thermoMixtureType::he,
        p,
        T
    );
}


template<class MixtureType, class BasicThermoType>
Foam::tmp<Foam::volScalarField::Internal>
Foam::BasicThermo<MixtureType, BasicThermoType>::he
(
    const volScalarField::Internal& p,
    const volScalarField::Internal& T
) const
{
    return volInternalScalarFieldProperty
    (
        "he",
        dimEnergy/dimMass,
        &MixtureType::thermoMixture,
        &MixtureType::thermoMixtureType::he,
        p,
        T
    );
}


template<class MixtureType, class BasicThermoType>
Foam::tmp<Foam::scalarField>
Foam::BasicThermo<MixtureType, BasicThermoType>::he
(
    const scalarField& T,
    const labelList& cells
) const
{
    return cellSetProperty
    (
        &MixtureType::thermoMixture,
        &MixtureType::thermoMixtureType::he,
        cells,
        cellSetScalarList(this->p_, cells),
        T
    );
}


template<class MixtureType, class BasicThermoType>
Foam::tmp<Foam::scalarField>
Foam::BasicThermo<MixtureType, BasicThermoType>::he
(
    const scalarField& T,
    const label patchi
) const
{
    return patchFieldProperty
    (
        &MixtureType::thermoMixture,
        &MixtureType::thermoMixtureType::he,
        patchi,
        this->p_.boundaryField()[patchi],
        T
    );
}


template<class MixtureType, class BasicThermoType>
Foam::tmp<Foam::scalarField>
Foam::BasicThermo<MixtureType, BasicThermoType>::he
(
    const scalarField& T,
    const fvSource& source
) const
{
    return fieldSourceProperty
    (
        &MixtureType::thermoMixture,
        &MixtureType::thermoMixtureType::he,
        source,
        cellSetScalarList(this->p_, source.cells()),
        T
    );
}


template<class MixtureType, class BasicThermoType>
Foam::tmp<Foam::volScalarField>
Foam::BasicThermo<MixtureType, BasicThermoType>::hs() const
{
    return volScalarFieldProperty
    (
        "hs",
        dimEnergy/dimMass,
        &MixtureType::thermoMixture,
        &MixtureType::thermoMixtureType::hs,
        this->p_,
        this->T_
    );
}


template<class MixtureType, class BasicThermoType>
Foam::tmp<Foam::volScalarField>
Foam::BasicThermo<MixtureType, BasicThermoType>::hs
(
    const volScalarField& p,
    const volScalarField& T
) const
{
    return volScalarFieldProperty
    (
        "hs",
        dimEnergy/dimMass,
        &MixtureType::thermoMixture,
        &MixtureType::thermoMixtureType::hs,
        p,
        T
    );
}


template<class MixtureType, class BasicThermoType>
Foam::tmp<Foam::volScalarField::Internal>
Foam::BasicThermo<MixtureType, BasicThermoType>::hs
(
    const volScalarField::Internal& p,
    const volScalarField::Internal& T
) const
{
    return volInternalScalarFieldProperty
    (
        "hs",
        dimEnergy/dimMass,
        &MixtureType::thermoMixture,
        &MixtureType::thermoMixtureType::hs,
        p,
        T
    );
}


template<class MixtureType, class BasicThermoType>
Foam::tmp<Foam::scalarField>
Foam::BasicThermo<MixtureType, BasicThermoType>::hs
(
    const scalarField& T,
    const labelList& cells
) const
{
    return cellSetProperty
    (
        &MixtureType::thermoMixture,
        &MixtureType::thermoMixtureType::hs,
        cells,
        cellSetScalarList(this->p_, cells),
        T
    );
}


template<class MixtureType, class BasicThermoType>
Foam::tmp<Foam::scalarField>
Foam::BasicThermo<MixtureType, BasicThermoType>::hs
(
    const scalarField& T,
    const label patchi
) const
{
    return patchFieldProperty
    (
        &MixtureType::thermoMixture,
        &MixtureType::thermoMixtureType::hs,
        patchi,
        this->p_.boundaryField()[patchi],
        T
    );
}


template<class MixtureType, class BasicThermoType>
Foam::tmp<Foam::volScalarField>
Foam::BasicThermo<MixtureType, BasicThermoType>::ha() const
{
    return volScalarFieldProperty
    (
        "ha",
        dimEnergy/dimMass,
        &MixtureType::thermoMixture,
        &MixtureType::thermoMixtureType::ha,
        this->p_,
        this->T_
    );
}


template<class MixtureType, class BasicThermoType>
Foam::tmp<Foam::volScalarField>
Foam::BasicThermo<MixtureType, BasicThermoType>::ha
(
    const volScalarField& p,
    const volScalarField& T
) const
{
    return volScalarFieldProperty
    (
        "ha",
        dimEnergy/dimMass,
        &MixtureType::thermoMixture,
        &MixtureType::thermoMixtureType::ha,
        p,
        T
    );
}


template<class MixtureType, class BasicThermoType>
Foam::tmp<Foam::volScalarField::Internal>
Foam::BasicThermo<MixtureType, BasicThermoType>::ha
(
    const volScalarField::Internal& p,
    const volScalarField::Internal& T
) const
{
    return volInternalScalarFieldProperty
    (
        "ha",
        dimEnergy/dimMass,
        &MixtureType::thermoMixture,
        &MixtureType::thermoMixtureType::ha,
        p,
        T
    );
}


template<class MixtureType, class BasicThermoType>
Foam::tmp<Foam::scalarField>
Foam::BasicThermo<MixtureType, BasicThermoType>::ha
(
    const scalarField& T,
    const labelList& cells
) const
{
    return cellSetProperty
    (
        &MixtureType::thermoMixture,
        &MixtureType::thermoMixtureType::ha,
        cells,
        cellSetScalarList(this->p_, cells),
        T
    );
}


template<class MixtureType, class BasicThermoType>
Foam::tmp<Foam::scalarField>
Foam::BasicThermo<MixtureType, BasicThermoType>::ha
(
    const scalarField& T,
    const label patchi
) const
{
    return patchFieldProperty
    (
        &MixtureType::thermoMixture,
        &MixtureType::thermoMixtureType::ha,
        patchi,
        this->p_.boundaryField()[patchi],
        T
    );
}


template<class MixtureType, class BasicThermoType>
Foam::tmp<Foam::scalarField>
Foam::BasicThermo<MixtureType, BasicThermoType>::Cp
(
    const scalarField& T,
    const label patchi
) const
{
    return patchFieldProperty
    (
        &MixtureType::thermoMixture,
        &MixtureType::thermoMixtureType::Cp,
        patchi,
        this->p_.boundaryField()[patchi],
        T
    );
}


template<class MixtureType, class BasicThermoType>
Foam::tmp<Foam::scalarField>
Foam::BasicThermo<MixtureType, BasicThermoType>::Cv
(
    const scalarField& T,
    const label patchi
) const
{
    return patchFieldProperty
    (
        &MixtureType::thermoMixture,
        &MixtureType::thermoMixtureType::Cv,
        patchi,
        this->p_.boundaryField()[patchi],
        T
    );
}


template<class MixtureType, class BasicThermoType>
Foam::tmp<Foam::scalarField>
Foam::BasicThermo<MixtureType, BasicThermoType>::Cpv
(
    const scalarField& T,
    const label patchi
) const
{
    if (MixtureType::thermoType::enthalpy())
    {
        return Cp(T, patchi);
    }
    else
    {
        return Cv(T, patchi);
    }
}


template<class MixtureType, class BasicThermoType>
Foam::tmp<Foam::volScalarField>
Foam::BasicThermo<MixtureType, BasicThermoType>::The
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
        &MixtureType::thermoMixture,
        &MixtureType::thermoMixtureType::The,
        h,
        p,
        T0
    );
}


template<class MixtureType, class BasicThermoType>
Foam::tmp<Foam::scalarField>
Foam::BasicThermo<MixtureType, BasicThermoType>::The
(
    const scalarField& h,
    const scalarField& T0,
    const labelList& cells
) const
{
    return cellSetProperty
    (
        &MixtureType::thermoMixture,
        &MixtureType::thermoMixtureType::The,
        cells,
        h,
        cellSetScalarList(this->p_, cells),
        T0
    );
}


template<class MixtureType, class BasicThermoType>
Foam::tmp<Foam::scalarField>
Foam::BasicThermo<MixtureType, BasicThermoType>::The
(
    const scalarField& h,
    const scalarField& T0,
    const label patchi
) const
{
    return patchFieldProperty
    (
        &MixtureType::thermoMixture,
        &MixtureType::thermoMixtureType::The,
        patchi,
        h,
        this->p_.boundaryField()[patchi],
        T0
    );
}


template<class MixtureType, class BasicThermoType>
bool Foam::BasicThermo<MixtureType, BasicThermoType>::read()
{
    if (physicalProperties::read())
    {
        MixtureType::read(*this);
        BasicThermoType::read(*this);
        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
