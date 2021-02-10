/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2021 OpenFOAM Foundation
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
            cellSetOption,
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
    volumeMode_ = volumeModeNames_.read(coeffs_.lookup("volumeMode"));

    // Set volume normalisation
    switch (volumeMode_)
    {
        case volumeMode::absolute:
            VDash_ = V();
            break;
        case volumeMode::specific:
            VDash_ = 1;
            break;
    }

    // Set field source terms
    fieldSp_.clear();
    fieldSu_.clear();
    forAllConstIter(dictionary, coeffs_.subDict("sources"), iter)
    {
        fieldSu_.set
        (
            iter().keyword(),
            objectFunction1::New<VolField>
            (
                "explicit",
                iter().dict(),
                iter().keyword(),
                mesh_,
                false
            ).ptr()
        );
        fieldSp_.set
        (
            iter().keyword(),
            Function1<scalar>::New
            (
                "implicit",
                iter().dict()
            ).ptr()
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
            << ">::addSup for source " << name_ << endl;
    }

    const scalar t = mesh_.time().value();

    const GeometricField<Type, fvPatchField, volMesh>& psi = eqn.psi();

    typename GeometricField<Type, fvPatchField, volMesh>::Internal Su
    (
        IOobject
        (
            name_ + fieldName + "Su",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensioned<Type>
        (
            "zero",
            eqn.dimensions()/dimVolume,
            Zero
        ),
        false
    );

    // Explicit source function for the field
    UIndirectList<Type>(Su, cells()) =
        fieldSu_[fieldName]->value<Type>(t)/VDash_;

    volScalarField::Internal Sp
    (
        IOobject
        (
            name_ + fieldName + "Sp",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensioned<scalar>
        (
            "zero",
            Su.dimensions()/psi.dimensions(),
            0
        ),
        false
    );

    // Implicit source function for the field
    UIndirectList<scalar>(Sp, cells()) =
        fieldSp_[fieldName]->value(t)/VDash_;

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
    cellSetOption(name, modelType, dict, mesh),
    volumeMode_(volumeMode::absolute),
    VDash_(1)
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


FOR_ALL_FIELD_TYPES(IMPLEMENT_FV_OPTION_ADD_SUP, semiImplicitSource);


FOR_ALL_FIELD_TYPES(IMPLEMENT_FV_OPTION_ADD_RHO_SUP, semiImplicitSource);


FOR_ALL_FIELD_TYPES(IMPLEMENT_FV_OPTION_ADD_ALPHA_RHO_SUP, semiImplicitSource);


bool Foam::fv::semiImplicitSource::read(const dictionary& dict)
{
    if (cellSetOption::read(dict))
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
