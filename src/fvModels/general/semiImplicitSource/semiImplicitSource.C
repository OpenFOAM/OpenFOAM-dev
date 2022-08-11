/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2022 OpenFOAM Foundation
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

void Foam::fv::semiImplicitSource::readCoeffs()
{
    // Get the volume mode
    volumeMode_ = volumeModeNames_.read(coeffs().lookup("volumeMode"));

    // Set field source terms
    fieldSu_.clear();
    fieldSp_.clear();
    forAllConstIter(dictionary, coeffs().subDict("sources"), iter)
    {
        fieldSu_.set
        (
            iter().keyword(),
            new unknownTypeFunction1("explicit", iter().dict())
        );
        fieldSp_.set
        (
            iter().keyword(),
            Function1<scalar>::New("implicit", iter().dict()).ptr()
        );
    }
}


template<class Type>
void Foam::fv::semiImplicitSource::addSupType
(
    fvMatrix<Type>& eqn,
    const word& fieldName
) const
{
    if (debug)
    {
        Info<< "semiImplicitSource<" << pTraits<Type>::typeName
            << ">::addSup for source " << name() << endl;
    }

    const scalar t = mesh().time().userTimeValue();

    const GeometricField<Type, fvPatchField, volMesh>& psi = eqn.psi();

    typename GeometricField<Type, fvPatchField, volMesh>::Internal Su
    (
        IOobject
        (
            name() + fieldName + "Su",
            mesh().time().timeName(),
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
        fieldSu_[fieldName]->value<Type>(t)/VDash;

    volScalarField::Internal Sp
    (
        IOobject
        (
            name() + fieldName + "Sp",
            mesh().time().timeName(),
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
        fieldSp_[fieldName]->value(t)/VDash;

    eqn += Su + fvm::SuSp(Sp, psi);
}


template<class Type>
void Foam::fv::semiImplicitSource::addSupType
(
    const volScalarField& rho,
    fvMatrix<Type>& eqn,
    const word& fieldName
) const
{
    return this->addSup(eqn, fieldName);
}


template<class Type>
void Foam::fv::semiImplicitSource::addSupType
(
    const volScalarField& alpha,
    const volScalarField& rho,
    fvMatrix<Type>& eqn,
    const word& fieldName
) const
{
    return this->addSup(eqn, fieldName);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::semiImplicitSource::semiImplicitSource
(
    const word& name,
    const word& modelType,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    fvModel(name, modelType, dict, mesh),
    set_(mesh, coeffs()),
    volumeMode_(volumeMode::absolute)
{
    readCoeffs();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fv::semiImplicitSource::~semiImplicitSource()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::wordList Foam::fv::semiImplicitSource::addSupFields() const
{
    return fieldSu_.toc();
}


FOR_ALL_FIELD_TYPES(IMPLEMENT_FV_MODEL_ADD_SUP, fv::semiImplicitSource);


FOR_ALL_FIELD_TYPES(IMPLEMENT_FV_MODEL_ADD_RHO_SUP, fv::semiImplicitSource);


FOR_ALL_FIELD_TYPES
(
    IMPLEMENT_FV_MODEL_ADD_ALPHA_RHO_SUP,
    fv::semiImplicitSource
);


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
        set_.read(coeffs());
        readCoeffs();
        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
