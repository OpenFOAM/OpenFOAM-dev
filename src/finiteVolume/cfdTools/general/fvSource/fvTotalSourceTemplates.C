/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2021-2025 OpenFOAM Foundation
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

#include "fvTotalSource.H"
#include "volFields.H"

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class Type>
void Foam::fvTotalSource::addSupType
(
    const VolField<Type>& field,
    fvMatrix<Type>& eqn
) const
{
    DebugInFunction
        << "field=" << field.name()
        << ", eqnField=" << eqn.psi().name() << endl;

    const labelUList cells = this->cells();

    const scalarField S
    (
        scalarField(mesh().V(), cells)/V()*this->S().value()
    );

    // Check the dimensions
    eqn.dimensions() = this->S().dimensions()*field.dimensions();

    if (&field == &eqn.psi())
    {
        // Get the field source coefficients
        const Field<Type> sourceCoeff
        (
            field.sources()[name()].sourceCoeff(*this, S, cells)
        );
        const scalarField internalCoeff
        (
            field.sources()[name()].internalCoeff(*this, S, cells)
        );

        // Apply the source
        Field<Type>& eqnSource = eqn.source();
        scalarField& eqnDiag = eqn.diag();
        forAll(cells, i)
        {
            eqnSource[cells[i]] -= S[i]*sourceCoeff[i];
            eqnDiag[cells[i]] += S[i]*internalCoeff[i];
        }
    }
    else
    {
        // Get the field source value
        const Field<Type> value
        (
            field.sources()[name()].value(*this, S, cells)
        );

        // Apply the source
        Field<Type>& eqnSource = eqn.source();
        forAll(cells, i)
        {
            eqnSource[cells[i]] -= S[i]*value[i];
        }
    }
}


template<class Type>
void Foam::fvTotalSource::addSupType
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
void Foam::fvTotalSource::addSupType
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
