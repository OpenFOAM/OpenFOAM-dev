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

#include "massTransferBase.H"
#include "fvMatrices.H"

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class Type>
const Foam::VolField<Type>& Foam::fv::massTransferBase::otherField
(
    const VolField<Type>& field
) const
{
    const label i = index(phaseNames(), field.group());

    const word otherFieldName =
        IOobject::groupName(field.member(), phaseNames()[!i]);

    if (mesh().foundObject<VolField<Type>>(otherFieldName))
    {
        return mesh().lookupObject<VolField<Type>>(otherFieldName);
    }
    else
    {
        return NullObjectRef<VolField<Type>>();
    }
}


template<class Type>
void Foam::fv::massTransferBase::addSupType
(
    const VolField<Type>& field,
    fvMatrix<Type>& eqn
) const
{
    DebugInFunction
        << "field=" << field.name()
        << ", eqnField=" << eqn.psi().name() << endl;

    // Volume-weighted mixture property equation
    if (field.group() == word::null)
    {
        eqn -= fvm::SuSp((1/rho(0) - 1/rho(1))*mDot(), eqn.psi());
    }
    // Not recognised. Fail.
    else
    {
        FatalErrorInFunction
            << "Cannot add a phase transfer for field " << field.name()
            << exit(FatalError);
    }
}


template<class Type>
void Foam::fv::massTransferBase::addSupType
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

    const label i = index(alphaNames(), alphaOrRho.name());

    // Incompressible property equation
    if
    (
        i != -1
     && i == index(phaseNames(), field.group())
    )
    {
        // Try and lookup the corresponding field in the other phase
        const VolField<Type>& otherField = this->otherField(field);

        // If a corresponding field exists then transfer its value
        if (!isNull(otherField))
        {
            const volScalarField::Internal S(this->S(field.name())/rho(i));

            eqn += posPart(S)*otherField;

            if (&field == &eqn.psi())
            {
                eqn += fvm::Sp(negPart(S), eqn.psi());
            }
            else
            {
                eqn += negPart(S)*field;
            }
        }
        // Otherwise, get the base class to create the source. This will
        // require a field source specification.
        else
        {
            fvMatrix<Type> rhoEqn(eqn.psi(), dimDensity*eqn.dimensions());

            fvSpecificSource::addSupType(alphaOrRho, field, rhoEqn);

            eqn += rhoEqn/rho(i);
        }
    }
    // Mass-weighted mixture property equation
    else if
    (
        alphaOrRho.group() == word::null
     && alphaOrRho.dimensions() == dimDensity
     && field.group() == word::null
    )
    {
        // Nothing to do. Standard transport terms already represent the
        // transfer. No source term is needed.
    }
    // Not recognised. Fail.
    else
    {
        FatalErrorInFunction
            << "Cannot add a phase transfer for field " << field.name()
            << exit(FatalError);
    }
}


template<class Type>
void Foam::fv::massTransferBase::addSupType
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

    const label i = index(alphaNames(), alpha.name());

    // Compressible property equation
    if
    (
        i != -1
     && i == index(rhoNames(), rho.name())
     && i == index(phaseNames(), field.group())
    )
    {
        // Try and lookup the corresponding field in the other phase
        const VolField<Type>& otherField = this->otherField(field);

        // If a corresponding field exists then transfer its value
        if (!isNull(otherField))
        {
            const volScalarField::Internal S(this->S(field.name()));

            eqn += posPart(S)*otherField;

            if (&field == &eqn.psi())
            {
                eqn += fvm::Sp(negPart(S), eqn.psi());
            }
            else
            {
                eqn += negPart(S)*field;
            }
        }
        // Otherwise, get the base class to create the source. This will
        // require a field source specification.
        else
        {
            fvSpecificSource::addSupType(alpha, rho, field, eqn);
        }
    }
    // Not recognised. Fail.
    else
    {
        FatalErrorInFunction
            << "Cannot add a phase transfer for field " << field.name()
            << exit(FatalError);
    }
}


// ************************************************************************* //
