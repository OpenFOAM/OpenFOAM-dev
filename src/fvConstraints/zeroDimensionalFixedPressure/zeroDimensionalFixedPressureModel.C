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
    const VolField<Type>& field,
    fvMatrix<Type>& eqn
) const
{
    FatalErrorInFunction
        << "Cannot add a fixed pressure source for field " << field.name()
        << " to equation for " << eqn.psi().name() << " because this field's "
        << "equation was not recognised as being in mass-conservative form"
        << exit(FatalError);
}


void Foam::fv::zeroDimensionalFixedPressureModel::addSupType
(
    const volScalarField& rho,
    fvMatrix<scalar>& eqn
) const
{
    if (IOobject::member(rho.name()) == constraint().rhoName())
    {
        if (IOobject::member(eqn.psi().name()) == constraint().pName())
        {
            eqn += constraint().pEqnSource(rho, eqn);
        }
        else
        {
            eqn += constraint().massSource(rho());
        }
    }
    else
    {
        // This is actually an incompressible single-phase equation. Rho is
        // actually a property field. Fall back (and error).
        addSupType<scalar>(rho, eqn);
    }
}


template<class Type>
void Foam::fv::zeroDimensionalFixedPressureModel::addSupType
(
    const volScalarField& rho,
    const VolField<Type>& field,
    fvMatrix<Type>& eqn
) const
{
    if (&field != &eqn.psi())
    {
        FatalErrorInFunction
            << "Cannot add a fixed pressure source of field " << field.name()
            << " into an equation for field " << eqn.psi().name()
            << exit(FatalError);
    }

    eqn -= fvm::SuSp(-constraint().massSource(rho()), eqn.psi());
}


void Foam::fv::zeroDimensionalFixedPressureModel::addSupType
(
    const volScalarField& alpha,
    const volScalarField& rho,
    fvMatrix<scalar>& eqn
) const
{
    if (IOobject::member(rho.name()) == constraint().rhoName())
    {
        if (IOobject::member(eqn.psi().name()) == constraint().pName())
        {
            eqn += alpha()*constraint().pEqnSource(rho, eqn);
        }
        else
        {
            eqn += constraint().massSource(alpha(), rho());
        }
    }
    else
    {
        // This is actually a compressible single-phase equation. Alpha is
        // actually rho, and rho is actually a property field. Fall back.
        addSupType<scalar>(alpha, rho, eqn);
    }
}


template<class Type>
void Foam::fv::zeroDimensionalFixedPressureModel::addSupType
(
    const volScalarField& alpha,
    const volScalarField& rho,
    const VolField<Type>& field,
    fvMatrix<Type>& eqn
) const
{
    if (&field != &eqn.psi())
    {
        FatalErrorInFunction
            << "Cannot add a fixed pressure source of field " << field.name()
            << " into an equation for field " << eqn.psi().name()
            << exit(FatalError);
    }

    eqn -= fvm::SuSp(-constraint().massSource(alpha(), rho()), eqn.psi());
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
    IMPLEMENT_FV_MODEL_ADD_FIELD_SUP,
    fv::zeroDimensionalFixedPressureModel
)


FOR_ALL_FIELD_TYPES
(
    IMPLEMENT_FV_MODEL_ADD_RHO_FIELD_SUP,
    fv::zeroDimensionalFixedPressureModel
)


FOR_ALL_FIELD_TYPES
(
    IMPLEMENT_FV_MODEL_ADD_ALPHA_RHO_FIELD_SUP,
    fv::zeroDimensionalFixedPressureModel
)


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
