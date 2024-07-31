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

#include "semiImplicitSource.H"
#include "fvMesh.H"
#include "fvMatrices.H"
#include "fvmSup.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace fv
    {
        defineTypeNameAndDebug(semiImplicitSource, 0);

        addToRunTimeSelectionTable
        (
            fvModel,
            semiImplicitSource,
            dictionary
        );
    }

    template<>
    const char* NamedEnum<fv::semiImplicitSource::volumeMode, 2>::names[] =
    {
        "absolute",
        "specific"
    };
}

const Foam::NamedEnum<Foam::fv::semiImplicitSource::volumeMode, 2>
    Foam::fv::semiImplicitSource::volumeModeNames_;


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::fv::semiImplicitSource::readCoeffs(const dictionary& dict)
{
    // Get the volume mode
    volumeMode_ = volumeModeNames_.read(dict.lookup("volumeMode"));

    // Set field source terms
    fieldSu_.clear();
    fieldSp_.clear();
    forAllConstIter(dictionary, dict.subDict("sources"), iter)
    {
        fieldSu_.set
        (
            iter().keyword(),
            new unknownTypeFunction1
            (
                "explicit",
                mesh().time().userUnits(),
                iter().dict()
            )
        );
        fieldSp_.set
        (
            iter().keyword(),
            new unknownTypeFunction1
            (
                "implicit",
                mesh().time().userUnits(),
                iter().dict()
            )
        );
    }
}


template<class Type>
void Foam::fv::semiImplicitSource::addSupType
(
    const VolField<Type>& field,
    fvMatrix<Type>& eqn
) const
{
    // Set the value units for the functions
    fieldSu_[field.name()]->template setValueUnits<Type>
    (
        eqn.dimensions()
    );
    fieldSp_[field.name()]->template setValueUnits<scalar>
    (
        eqn.dimensions()/eqn.psi().dimensions()
    );

    const scalar t = mesh().time().value();

    const VolField<Type>& psi = eqn.psi();

    VolInternalField<Type> Su
    (
        IOobject
        (
            name() + field.name() + "Su",
            mesh().time().name(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh(),
        dimensioned<Type>
        (
            "zero",
            eqn.dimensions()/dimVolume,
            Zero
        ),
        false
    );

    // Set volume normalisation
    scalar VDash = NaN;
    switch (volumeMode_)
    {
        case volumeMode::absolute:
            VDash = set_.V();
            break;
        case volumeMode::specific:
            VDash = 1;
            break;
    }

    // Explicit source function for the field
    UIndirectList<Type>(Su, set_.cells()) =
        fieldSu_[field.name()]->template value<Type>(t)/VDash;

    volScalarField::Internal Sp
    (
        IOobject
        (
            name() + field.name() + "Sp",
            mesh().time().name(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh(),
        dimensioned<scalar>
        (
            "zero",
            Su.dimensions()/psi.dimensions(),
            0
        ),
        false
    );

    // Implicit source function for the field
    UIndirectList<scalar>(Sp, set_.cells()) =
        fieldSp_[field.name()]->template value<scalar>(t)/VDash;

    eqn += Su - fvm::SuSp(-Sp, psi);
}


template<class Type>
void Foam::fv::semiImplicitSource::addSupType
(
    const volScalarField& rho,
    const VolField<Type>& field,
    fvMatrix<Type>& eqn
) const
{
    return addSup(field, eqn);
}


template<class Type>
void Foam::fv::semiImplicitSource::addSupType
(
    const volScalarField& alpha,
    const volScalarField& rho,
    const VolField<Type>& field,
    fvMatrix<Type>& eqn
) const
{
    return addSup(field, eqn);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::semiImplicitSource::semiImplicitSource
(
    const word& name,
    const word& modelType,
    const fvMesh& mesh,
    const dictionary& dict
)
:
    fvModel(name, modelType, mesh, dict),
    set_(mesh, coeffs(dict)),
    volumeMode_(volumeMode::absolute)
{
    readCoeffs(coeffs(dict));
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fv::semiImplicitSource::~semiImplicitSource()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::wordList Foam::fv::semiImplicitSource::addSupFields() const
{
    return fieldSu_.toc();
}


FOR_ALL_FIELD_TYPES
(
    IMPLEMENT_FV_MODEL_ADD_FIELD_SUP,
    fv::semiImplicitSource
)


FOR_ALL_FIELD_TYPES
(
    IMPLEMENT_FV_MODEL_ADD_RHO_FIELD_SUP,
    fv::semiImplicitSource
)


FOR_ALL_FIELD_TYPES
(
    IMPLEMENT_FV_MODEL_ADD_ALPHA_RHO_FIELD_SUP,
    fv::semiImplicitSource
)


bool Foam::fv::semiImplicitSource::movePoints()
{
    set_.movePoints();
    return true;
}


void Foam::fv::semiImplicitSource::topoChange(const polyTopoChangeMap& map)
{
    set_.topoChange(map);
}


void Foam::fv::semiImplicitSource::mapMesh(const polyMeshMap& map)
{
    set_.mapMesh(map);
}


void Foam::fv::semiImplicitSource::distribute
(
    const polyDistributionMap& map
)
{
    set_.distribute(map);
}


bool Foam::fv::semiImplicitSource::read(const dictionary& dict)
{
    if (fvModel::read(dict))
    {
        set_.read(coeffs(dict));
        readCoeffs(coeffs(dict));
        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
