/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2026 OpenFOAM Foundation
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

Typedef
    Foam::tensorField

Description
    Specialisation of Field\<T\> for tensor.

SourceFiles
    tensorField.C

\*---------------------------------------------------------------------------*/

#ifndef tensorField_H
#define tensorField_H

#include "vectorField.H"
#include "tensorExpressions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
inline auto Foam::Field<Type>::T() const
{
    return Foam::T(*this);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

typedef Field<tensor> tensorField;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void inv
(
    Field<tensor>& tf,
    const UList<tensor>& tf1,
    const Vector<label>& solutionD
);

tmp<tensorField> inv
(
    const Field<tensor>& tf,
    const Vector<label>& solutionD
);

tmp<tensorField> inv
(
    const SubField<tensor>& tf,
    const Vector<label>& solutionD
);

tmp<tensorField> inv
(
    const tmp<tensorField>& tf,
    const Vector<label>& solutionD
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
