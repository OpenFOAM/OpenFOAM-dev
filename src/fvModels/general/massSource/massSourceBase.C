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

#include "massSourceBase.H"
#include "fvMatrices.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(massSourceBase, 0);
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::fv::massSourceBase::readCoeffs()
{
    rhoName_ =
        coeffs().lookupOrDefault<word>
        (
            "rho",
            IOobject::groupName("rho", phaseName())
        );
}


template<class Type>
void Foam::fv::massSourceBase::addSupType
(
    const VolField<Type>& field,
    fvMatrix<Type>& eqn
) const
{
    DebugInFunction
        << "field=" << field.name()
        << ", eqnField=" << eqn.psi().name() << endl;

    FatalErrorInFunction
        << "Cannot add a mass source for field " << field.name()
        << " to equation for " << eqn.psi().name() << " because this field's "
        << "equation was not recognised as being in mass-conservative form"
        << exit(FatalError);
}


void Foam::fv::massSourceBase::addSupType
(
    const volScalarField& rhoOrField,
    fvMatrix<scalar>& eqn
) const
{
    DebugInFunction
        << "rhoOrField=" << rhoOrField.name()
        << ", eqnField=" << eqn.psi().name() << endl;

    // Single-phase continuity equation
    if (phaseName() == word::null && rhoOrField.name() == rhoName_)
    {
        fvTotalSource::addSource(eqn);
    }
    // Not recognised. Fail.
    else
    {
        addSupType<scalar>(rhoOrField, eqn);
    }
}


template<class Type>
void Foam::fv::massSourceBase::addSupType
(
    const volScalarField& rho,
    const VolField<Type>& field,
    fvMatrix<Type>& eqn
) const
{
    DebugInFunction
        << "rho=" << rho.name()
        << ", field=" << field.name()
        << ", eqnField=" << eqn.psi().name() << endl;

    // Single-phase property equation
    if (phaseName() == word::null && rho.name() == rhoName_)
    {
        fvTotalSource::addSupType(rho, field, eqn);
    }
    // Multiphase mass-weighted mixture property equation
    else if
    (
        phaseName() != word::null
     && rho.group() == word::null
     && rho.dimensions() == dimDensity
     && field.group() == word::null
    )
    {
        fvTotalSource::addSupType(rho, field, eqn);
    }
    // Not recognised. Fail.
    else
    {
        addSupType<Type>(field, eqn);
    }
}


void Foam::fv::massSourceBase::addSupType
(
    const volScalarField& alphaOrRho,
    const volScalarField& rhoOrField,
    fvMatrix<scalar>& eqn
) const
{
    DebugInFunction
        << "alphaOrRho=" << alphaOrRho.name()
        << ", rhoOrField=" << rhoOrField.name()
        << ", eqnField=" << eqn.psi().name() << endl;

    // Multiphase continuity equation
    if (rhoOrField.name() == rhoName_)
    {
        fvTotalSource::addSource(eqn);
    }
    // Try the general type method
    else
    {
        addSupType<scalar>(alphaOrRho, rhoOrField, eqn);
    }
}


template<class Type>
void Foam::fv::massSourceBase::addSupType
(
    const volScalarField& alpha,
    const volScalarField& rho,
    const VolField<Type>& field,
    fvMatrix<Type>& eqn
) const
{
    DebugInFunction
        << "alpha=" << alpha.name()
        << ", rho=" << rho.name()
        << ", field=" << field.name()
        << ", eqnField=" << eqn.psi().name() << endl;

    // Multiphase property equation
    fvTotalSource::addSupType(alpha, rho, field, eqn);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::massSourceBase::massSourceBase
(
    const word& name,
    const word& modelType,
    const fvMesh& mesh,
    const dictionary& dict
)
:
    fvTotalSource(name, modelType, mesh, dict),
    rhoName_()
{
    readCoeffs();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fv::massSourceBase::addSup(fvMatrix<scalar>& eqn) const
{
    DebugInFunction
        << "eqnField=" << eqn.psi().name() << endl;

    FatalErrorInFunction
        << "Field-less mass sources are not possible"
        << exit(FatalError);
}


FOR_ALL_FIELD_TYPES(IMPLEMENT_FV_MODEL_ADD_FIELD_SUP, fv::massSourceBase);


FOR_ALL_FIELD_TYPES(IMPLEMENT_FV_MODEL_ADD_RHO_FIELD_SUP, fv::massSourceBase);


FOR_ALL_FIELD_TYPES
(
    IMPLEMENT_FV_MODEL_ADD_ALPHA_RHO_FIELD_SUP,
    fv::massSourceBase
);


bool Foam::fv::massSourceBase::read(const dictionary& dict)
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
