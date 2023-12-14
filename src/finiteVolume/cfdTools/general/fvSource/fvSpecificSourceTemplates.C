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

#include "fvSpecificSource.H"
#include "volFields.H"
#include "fvmSup.H"

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class Type>
void Foam::fvSpecificSource::addSupType
(
    const VolField<Type>& field,
    fvMatrix<Type>& eqn
) const
{
    DebugInFunction
        << "field=" << field.name()
        << ", eqnField=" << eqn.psi().name() << endl;

    tmp<volScalarField::Internal> tS(S(field.name()));

    if (&field == &eqn.psi())
    {
        // Get the field source coefficients
        tmp<typename VolField<Type>::Internal> sourceCoeff =
            VolField<Type>::Internal::New
            (
                "sourceCoeff",
                mesh(),
                field.dimensions(),
                field.sources()[name()].sourceCoeff(*this)
            );
        tmp<typename volScalarField::Internal> internalCoeff =
            volScalarField::Internal::New
            (
                "internalCoeff",
                mesh(),
                dimless,
                field.sources()[name()].internalCoeff(*this)
            );

        // Apply the source
        eqn += tS()*sourceCoeff + fvm::Sp(tS()*internalCoeff, eqn.psi());
    }
    else
    {
        // Get the field source value
        tmp<typename VolField<Type>::Internal> value =
            VolField<Type>::Internal::New
            (
                "value",
                mesh(),
                field.dimensions(),
                field.sources()[name()].value(*this)
            );

        // Apply the source
        eqn += tS*value;
    }
}


template<class Type>
void Foam::fvSpecificSource::addSupType
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

    addSupType(field, eqn);
}


template<class Type>
void Foam::fvSpecificSource::addSupType
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

    addSupType(rho, field, eqn);
}


// ************************************************************************* //
