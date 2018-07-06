/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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

InClass
    Foam::fvMatrix

Description
    A scalar instance of fvMatrix

SourceFiles
    fvScalarMatrix.C

\*---------------------------------------------------------------------------*/

#ifndef fvScalarMatrix_H
#define fvScalarMatrix_H

#include "fvMatrix.H"
#include "fvMatricesFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Set reference level for a component of the solution
// on a given patch face
template<>
void fvMatrix<scalar>::setComponentReference
(
    const label patchi,
    const label facei,
    const direction,
    const scalar value
);

template<>
autoPtr<fvMatrix<scalar>::fvSolver> fvMatrix<scalar>::solver
(
    const dictionary&
);

template<>
solverPerformance fvMatrix<scalar>::fvSolver::solve
(
    const dictionary&
);

template<>
solverPerformance fvMatrix<scalar>::solveSegregated
(
    const dictionary&
);

template<>
tmp<scalarField> fvMatrix<scalar>::residual() const;

template<>
tmp<volScalarField> fvMatrix<scalar>::H() const;

template<>
tmp<volScalarField> fvMatrix<scalar>::H1() const;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
