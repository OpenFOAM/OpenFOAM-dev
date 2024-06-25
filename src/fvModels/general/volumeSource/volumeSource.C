/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2021-2024 OpenFOAM Foundation
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

#include "volumeSource.H"
#include "fvMatrices.H"
#include "basicThermo.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(volumeSource, 0);
    addToRunTimeSelectionTable(fvModel, volumeSource, dictionary);
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::fv::volumeSource::readCoeffs()
{
    alphaName_ =
        phaseName() == word::null
      ? word::null
      : coeffs().lookupOrDefault<word>
        (
            "alpha",
            IOobject::groupName("alpha", phaseName())
        );

    setPtr_->read(coeffs());

    volumetricFlowRate_.reset
    (
        Function1<scalar>::New
        (
            "volumetricFlowRate",
            mesh().time().userUnits(),
            dimVolume/dimTime,
            coeffs()
        ).ptr()
    );
}


template<class Type>
void Foam::fv::volumeSource::addSupType
(
    const VolField<Type>& field,
    fvMatrix<Type>& eqn
) const
{
    DebugInFunction
        << "field=" << field.name()
        << ", eqnField=" << eqn.psi().name() << endl;

    // Single-phase property equation
    if (phaseName() == word::null && field.group() == word::null)
    {
        fvTotalSource::addSupType(field, eqn);
    }
    // Multiphase volume-weighted mixture property equation (e.g., a turbulence
    // equation if running standard incompressible transport modelling in the
    // incompressibleVoF solver)
    else if (phaseName() != word::null && field.group() == word::null)
    {
        fvTotalSource::addSupType(field, eqn);
    }
    // Not recognised. Fail.
    else
    {
        const volScalarField& null = NullObjectRef<volScalarField>();
        addSupType(null, null, field, eqn);
    }
}


void Foam::fv::volumeSource::addSupType
(
    const volScalarField& alphaOrField,
    fvMatrix<scalar>& eqn
) const
{
    DebugInFunction
        << "alphaOrField=" << alphaOrField.name()
        << ", eqnField=" << eqn.psi().name() << endl;

    // Multiphase continuity equation
    if (phaseName() != word::null && alphaOrField.name() == alphaName_)
    {
        fvTotalSource::addSource(eqn);
    }
    // Try the general type method
    else
    {
        addSupType<scalar>(alphaOrField, eqn);
    }
}


template<class Type>
void Foam::fv::volumeSource::addSupType
(
    const volScalarField& alphaOrRho,
    const VolField<Type>& field,
    fvMatrix<Type>& eqn
) const
{
    DebugInFunction
        << "alphaOrRho=" << alphaOrRho.name()
        << ", field=" << field.name()
        << ", eqnField=" << eqn.psi().name() << endl;

    // Multiphase property equation (e.g., a turbulence equation if running
    // two-phase transport modelling in the incompressibleVoF solver)
    if (phaseName() != word::null && alphaOrRho.name() == alphaName_)
    {
        fvTotalSource::addSupType(field, eqn);
    }
    // Multiphase mass-weighted mixture property equation (e.g., the momentum
    // equation in the incompressibleVoF solver)
    else if
    (
        phaseName() != word::null
     && alphaOrRho.group() == word::null
     && alphaOrRho.dimensions() == dimDensity
     && field.group() == word::null
    )
    {
        // First we construct the volumetric source...
        fvMatrix<Type> volEqn(eqn.psi(), eqn.dimensions()/dimDensity);
        fvTotalSource::addSupType(field, volEqn);

        // Then, to apply it to the mixture equation, we need to multiply by
        // the density of the phase of which this is a source. There is no
        // solver-agnostic interface, at present, that lets us obtain this
        // density. So, we read it from the physical properties file. This is
        // clunky, but it should work in all circumstances. This is what the
        // clouds fvModel does,
        const dimensionedScalar rhoi
        (
            "rho",
            dimDensity,
            mesh().lookupObject<IOdictionary>
            (
                IOobject::groupName
                (
                    physicalProperties::typeName,
                    phaseName()
                )
            )
        );
        eqn += rhoi*volEqn;
    }
    // Not recognised. Fail.
    else
    {
        const volScalarField& null = NullObjectRef<volScalarField>();
        addSupType(null, alphaOrRho, field, eqn);
    }
}


template<class Type>
void Foam::fv::volumeSource::addSupType
(
    const volScalarField& alpha,
    const volScalarField& rho,
    const VolField<Type>& field,
    fvMatrix<Type>& eqn
) const
{
    DebugInFunction
        << "alpha=" << (isNull(alpha) ? word::null : alpha.name())
        << ", rho=" << (isNull(rho) ? word::null : rho.name())
        << ", field=" << field.name()
        << ", eqnField=" << eqn.psi().name() << endl;

    FatalErrorInFunction
        << "Cannot add a volume source for field " << field.name()
        << " to equation for " << eqn.psi().name() << " because this field's "
        << "equation was not recognised as being in volume-conservative form"
        << exit(FatalError);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::volumeSource::volumeSource
(
    const word& name,
    const word& modelType,
    const fvMesh& mesh,
    const dictionary& dict
)
:
    fvTotalSource(name, modelType, mesh, dict),
    alphaName_(),
    setPtr_(new fvCellSet(mesh)),
    volumetricFlowRate_()
{
    readCoeffs();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::labelUList Foam::fv::volumeSource::cells() const
{
    return setPtr_->cells();
}


Foam::label Foam::fv::volumeSource::nCells() const
{
    return setPtr_->nCells();
}


Foam::scalar Foam::fv::volumeSource::V() const
{
    return setPtr_->V();
}


Foam::dimensionedScalar Foam::fv::volumeSource::S() const
{
    return
        dimensionedScalar
        (
            dimVolume/dimTime,
            volumetricFlowRate_->value(mesh().time().value())
        );
}


void Foam::fv::volumeSource::addSup(fvMatrix<scalar>& eqn) const
{
    DebugInFunction
        << "eqnField=" << eqn.psi().name() << endl;

    // Single-phase continuity equation
    fvTotalSource::addSource(eqn);
}


FOR_ALL_FIELD_TYPES(IMPLEMENT_FV_MODEL_ADD_FIELD_SUP, fv::volumeSource)


FOR_ALL_FIELD_TYPES(IMPLEMENT_FV_MODEL_ADD_RHO_FIELD_SUP, fv::volumeSource)


FOR_ALL_FIELD_TYPES
(
    IMPLEMENT_FV_MODEL_ADD_ALPHA_RHO_FIELD_SUP,
    fv::volumeSource
)


bool Foam::fv::volumeSource::movePoints()
{
    setPtr_->movePoints();
    return true;
}


void Foam::fv::volumeSource::topoChange(const polyTopoChangeMap& map)
{
    setPtr_->topoChange(map);
}


void Foam::fv::volumeSource::mapMesh(const polyMeshMap& map)
{
    setPtr_->mapMesh(map);
}


void Foam::fv::volumeSource::distribute(const polyDistributionMap& map)
{
    setPtr_->distribute(map);
}


bool Foam::fv::volumeSource::read(const dictionary& dict)
{
    if (fvTotalSource::read(dict))
    {
        readCoeffs();
        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
