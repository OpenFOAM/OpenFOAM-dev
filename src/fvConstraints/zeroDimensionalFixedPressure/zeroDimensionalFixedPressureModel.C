/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2023 OpenFOAM Foundation
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

#include "zeroDimensionalFixedPressureModel.H"
#include "zeroDimensionalFixedPressureConstraint.H"
#include "fvConstraints.H"
#include "fvmSup.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(zeroDimensionalFixedPressureModel, 0);
    addToRunTimeSelectionTable
    (
        fvModel,
        zeroDimensionalFixedPressureModel,
        dictionary
    );
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

const Foam::fv::zeroDimensionalFixedPressureConstraint&
Foam::fv::zeroDimensionalFixedPressureModel::constraint() const
{
    const fvConstraints& constraints = fvConstraints::New(mesh());

    forAll(constraints, i)
    {
        if (isA<zeroDimensionalFixedPressureConstraint>(constraints[i]))
        {
            return refCast<const zeroDimensionalFixedPressureConstraint>
            (
                constraints[i]
            );
        }
    }

    FatalErrorInFunction
        << "The " << typeName << " fvModel requires a corresponding "
        << zeroDimensionalFixedPressureConstraint::typeName << " fvConstraint"
        << exit(FatalError);

    return NullObjectRef<zeroDimensionalFixedPressureConstraint>();
}


template<class Type>
void Foam::fv::zeroDimensionalFixedPressureModel::addSupType
(
    fvMatrix<Type>& eqn,
    const word& fieldName
) const
{
    FatalErrorInFunction
        << "Cannot add a fixed pressure source to field " << fieldName
        << " because this field's equation is not in mass-conservative form"
        << exit(FatalError);
}


void Foam::fv::zeroDimensionalFixedPressureModel::addSupType
(
    fvMatrix<scalar>& eqn,
    const word& fieldName
) const
{
    if (IOobject::member(fieldName) == constraint().rhoName())
    {
        eqn += constraint().massSource(eqn.psi()());
    }
    else
    {
        addSupType<scalar>(eqn, fieldName); // error above
    }
}


template<class Type>
void Foam::fv::zeroDimensionalFixedPressureModel::addSupType
(
    const volScalarField& rho,
    fvMatrix<Type>& eqn,
    const word& fieldName
) const
{
    eqn -= fvm::SuSp(-constraint().massSource(rho()), eqn.psi());
}


void Foam::fv::zeroDimensionalFixedPressureModel::addSupType
(
    const volScalarField& rho,
    fvMatrix<scalar>& eqn,
    const word& fieldName
) const
{
    if (IOobject::member(fieldName) == constraint().rhoName())
    {
        if (IOobject::member(eqn.psi().name()) == constraint().pName())
        {
            eqn += constraint().pEqnSource(eqn);
        }
        else if (IOobject::member(eqn.psi().name()) == constraint().rhoName())
        {
            // Phase density equation. Argument names are misleading.
            const volScalarField& alpha = rho;
            const volScalarField& rho = eqn.psi();

            eqn += constraint().massSource(alpha(), rho());
        }
        else
        {
            FatalErrorInFunction
                << "Cannot add source for density field " << fieldName
                << " into an equation for " << eqn.psi().name()
                << exit(FatalError);
        }
    }
    else
    {
        addSupType<scalar>(rho, eqn, fieldName);
    }
}


template<class Type>
void Foam::fv::zeroDimensionalFixedPressureModel::addSupType
(
    const volScalarField& alpha,
    const volScalarField& rho,
    fvMatrix<Type>& eqn,
    const word& fieldName
) const
{
    eqn -= fvm::SuSp(-constraint().massSource(alpha(), rho()), eqn.psi());
}


void Foam::fv::zeroDimensionalFixedPressureModel::addSupType
(
    const volScalarField& alpha,
    const volScalarField& rho,
    fvMatrix<scalar>& eqn,
    const word& fieldName
) const
{
    if (IOobject::member(fieldName) == constraint().rhoName())
    {
        if (IOobject::member(eqn.psi().name()) == constraint().pName())
        {
            eqn += alpha*constraint().pEqnSource(eqn);
        }
        else if (IOobject::member(eqn.psi().name()) == constraint().rhoName())
        {
            FatalErrorInFunction
                << "Cannot add source for density field " << fieldName
                << " into a phase-conservative equation for "
                << eqn.psi().name() << exit(FatalError);
        }
        else
        {
            FatalErrorInFunction
                << "Cannot add source for density field " << fieldName
                << " into an equation for " << eqn.psi().name()
                << exit(FatalError);
        }
    }
    else
    {
        addSupType<scalar>(alpha, rho, eqn, fieldName);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::zeroDimensionalFixedPressureModel::zeroDimensionalFixedPressureModel
(
    const word& name,
    const word& modelType,
    const fvMesh& mesh,
    const dictionary& dict
)
:
    fvModel(name, modelType, mesh, dict)
{
    if (mesh.nGeometricD() != 0)
    {
        FatalIOErrorInFunction(dict)
            << "Zero-dimensional fvModel applied to a "
            << mesh.nGeometricD() << "-dimensional mesh"
            << exit(FatalIOError);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fv::zeroDimensionalFixedPressureModel::
~zeroDimensionalFixedPressureModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::fv::zeroDimensionalFixedPressureModel::addsSupToField
(
    const word& fieldName
) const
{
    return true;
}


FOR_ALL_FIELD_TYPES
(
    IMPLEMENT_FV_MODEL_ADD_SUP,
    fv::zeroDimensionalFixedPressureModel
);


FOR_ALL_FIELD_TYPES
(
    IMPLEMENT_FV_MODEL_ADD_RHO_SUP,
    fv::zeroDimensionalFixedPressureModel
);


FOR_ALL_FIELD_TYPES
(
    IMPLEMENT_FV_MODEL_ADD_ALPHA_RHO_SUP,
    fv::zeroDimensionalFixedPressureModel
);


bool Foam::fv::zeroDimensionalFixedPressureModel::movePoints()
{
    return true;
}


void Foam::fv::zeroDimensionalFixedPressureModel::topoChange
(
    const polyTopoChangeMap& map
)
{}


void Foam::fv::zeroDimensionalFixedPressureModel::mapMesh
(
    const polyMeshMap& map
)
{}


void Foam::fv::zeroDimensionalFixedPressureModel::distribute
(
    const polyDistributionMap& map
)
{}


// ************************************************************************* //
