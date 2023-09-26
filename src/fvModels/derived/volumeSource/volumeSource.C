/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2021-2023 OpenFOAM Foundation
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
    phaseName_ = coeffs().lookupOrDefault<word>("phase", word::null);

    alphaName_ =
        phaseName_ == word::null
      ? word::null
      : coeffs().lookupOrDefault<word>
        (
            "alpha",
            IOobject::groupName("alpha", phaseName_)
        );

    readSet();

    readFieldValues();

    volumetricFlowRate_.reset
    (
        Function1<scalar>::New("volumetricFlowRate", coeffs()).ptr()
    );
}


Foam::scalar Foam::fv::volumeSource::volumetricFlowRate() const
{
    return volumetricFlowRate_->value(mesh().time().userTimeValue());
}


template<class Type>
void Foam::fv::volumeSource::addSupType(fvMatrix<Type>& eqn) const
{
    FatalErrorInFunction
        << "Continuity sources for non-scalar types are not supported"
        << exit(FatalError);
}


void Foam::fv::volumeSource::addSupType(fvMatrix<scalar>& eqn) const
{
    const labelUList cells = set_.cells();

    const scalar volumetricFlowRate = this->volumetricFlowRate();

    // Continuity equation. Add the volumetric flow rate.
    forAll(cells, i)
    {
        eqn.source()[cells[i]] -=
            mesh().V()[cells[i]]/set_.V()*volumetricFlowRate;
    }
}


template<class Type>
void Foam::fv::volumeSource::addSupType
(
    const dimensionedScalar& oneOrRho,
    const VolField<Type>& field,
    fvMatrix<Type>& eqn
) const
{
    const labelUList cells = set_.cells();

    const scalar flowRate = this->volumetricFlowRate()*oneOrRho.value();

    // Property equation. If the source is positive, introduce the value
    // specified by the user. If negative, then sink the current internal value
    // using an implicit term.
    if (flowRate > 0)
    {
        const Type value =
            fieldValues_[field.name()]->template value<Type>
            (
                mesh().time().userTimeValue()
            );

        forAll(cells, i)
        {
            eqn.source()[cells[i]] -=
                mesh().V()[cells[i]]/set_.V()*flowRate*value;
        }
    }
    else
    {
        forAll(cells, i)
        {
            eqn.diag()[cells[i]] +=
                mesh().V()[cells[i]]/set_.V()*flowRate;
        }
    }
}


template<class Type>
void Foam::fv::volumeSource::addSupType
(
    const VolField<Type>& field,
    fvMatrix<Type>& eqn
) const
{
    // Property equation
    addSupType(dimensionedScalar(dimless, scalar(1)), field, eqn);
}


void Foam::fv::volumeSource::addSupType
(
    const volScalarField& field,
    fvMatrix<scalar>& eqn
) const
{
    // Multiphase continuity equation. Same source as single-phase case.
    if (field.name() == alphaName_)
    {
        addSupType(eqn);
        return;
    }

    // Property equation
    addSupType<scalar>(field, eqn);
}


template<class Type>
void Foam::fv::volumeSource::addSupType
(
    const volScalarField& rho,
    const VolField<Type>& field,
    fvMatrix<Type>& eqn
) const
{
    // Multiphase property equation (e.g., turbulence equation if running
    // two-phase transport modelling in the incompressibleVoF solver)
    if (rho.name() == alphaName_)
    {
        addSupType(dimensionedScalar(dimless, scalar(1)), field, eqn);
        return;
    }

    // Mixture property equation (e.g., the momentum equation in the
    // incompressibleVoF solver)...

    // We need to know the density of the phase of which this is a source in
    // order to create the relevant term. There is no solver-agnostic
    // interface, at present, that lets us do this. So, read the density from
    // the physical properties file. This is clunky, but it should work in all
    // circumstances. This is what the clouds fvModel does,
    const dimensionedScalar rhoi
    (
        "rho",
        dimDensity,
        mesh().lookupObject<IOdictionary>
        (
            IOobject::groupName
            (
                physicalProperties::typeName,
                phaseName_
            )
        )
    );

    addSupType(rhoi, field, eqn);
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
    FatalErrorInFunction
        << "Cannot add a volume source for field " << field.name()
        << " because this field's equation is not in volume-conservative form"
        << exit(FatalError);
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::fv::volumeSource::readSet()
{
    set_.read(coeffs());
}


void Foam::fv::volumeSource::readFieldValues()
{
    fieldValues_.clear();
    const dictionary& fieldCoeffs = coeffs().subDict("fieldValues");
    forAllConstIter(dictionary, fieldCoeffs, iter)
    {
        fieldValues_.set
        (
            iter().keyword(),
            new unknownTypeFunction1(iter().keyword(), fieldCoeffs)
        );
    }
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
    fvModel(name, modelType, mesh, dict),
    phaseName_(),
    set_(fvCellSet(mesh)),
    fieldValues_(),
    volumetricFlowRate_()
{
    readCoeffs();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::fv::volumeSource::addsSupToField(const word& fieldName) const
{
    const bool isMixture = IOobject::group(fieldName) == word::null;
    const bool isThisPhase = IOobject::group(fieldName) == phaseName_;

    if
    (
        (isMixture || isThisPhase)
     && volumetricFlowRate() > 0
     && !(fieldName == alphaName_)
     && !fieldValues_.found(fieldName)
    )
    {
        WarningInFunction
            << "No value supplied for field " << fieldName << " in "
            << type() << " fvModel " << name() << endl;

        return false;
    }

    return isMixture || isThisPhase;
}


Foam::wordList Foam::fv::volumeSource::addSupFields() const
{
    return fieldValues_.toc();
}


FOR_ALL_FIELD_TYPES
(
    IMPLEMENT_FV_MODEL_ADD_SUP,
    fv::volumeSource
)


FOR_ALL_FIELD_TYPES
(
    IMPLEMENT_FV_MODEL_ADD_FIELD_SUP,
    fv::volumeSource
)


FOR_ALL_FIELD_TYPES
(
    IMPLEMENT_FV_MODEL_ADD_RHO_FIELD_SUP,
    fv::volumeSource
)


FOR_ALL_FIELD_TYPES
(
    IMPLEMENT_FV_MODEL_ADD_ALPHA_RHO_FIELD_SUP,
    fv::volumeSource
)


bool Foam::fv::volumeSource::movePoints()
{
    set_.movePoints();
    return true;
}


void Foam::fv::volumeSource::topoChange(const polyTopoChangeMap& map)
{
    set_.topoChange(map);
}


void Foam::fv::volumeSource::mapMesh(const polyMeshMap& map)
{
    set_.mapMesh(map);
}


void Foam::fv::volumeSource::distribute(const polyDistributionMap& map)
{
    set_.distribute(map);
}


bool Foam::fv::volumeSource::read(const dictionary& dict)
{
    if (fvModel::read(dict))
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
