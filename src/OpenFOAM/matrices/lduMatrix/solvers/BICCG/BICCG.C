/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

#include "BICCG.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(BICCG, 0);

    lduMatrix::solver::addasymMatrixConstructorToTable<BICCG>
        addBICCGSymMatrixConstructorToTable_;
}


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

Foam::dictionary Foam::BICCG::solverDict
(
    const scalar tol,
    const scalar relTol
)
{
    dictionary dict(IStringStream("solver PBiCG; preconditioner DILU;")());
    dict.add("tolerance", tol);
    dict.add("relTol", relTol);

    return dict;
}


Foam::dictionary Foam::BICCG::solverDict
(
    Istream& is
)
{
    scalar tol(readScalar(is));
    scalar relTol(readScalar(is));

    return solverDict(tol, relTol);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::BICCG::BICCG
(
    const word& fieldName,
    const lduMatrix& matrix,
    const FieldField<Field, scalar>& interfaceBouCoeffs,
    const FieldField<Field, scalar>& interfaceIntCoeffs,
    const lduInterfaceFieldPtrsList& interfaces,
    const dictionary& solverControls
)
:
    PBiCG
    (
        fieldName,
        matrix,
        interfaceBouCoeffs,
        interfaceIntCoeffs,
        interfaces,
        solverControls
    )
{}


Foam::BICCG::BICCG
(
    const word& fieldName,
    const lduMatrix& matrix,
    const FieldField<Field, scalar>& interfaceBouCoeffs,
    const FieldField<Field, scalar>& interfaceIntCoeffs,
    const lduInterfaceFieldPtrsList& interfaces,
    const scalar tolerance,
    const scalar relTol
)
:
    PBiCG
    (
        fieldName,
        matrix,
        interfaceBouCoeffs,
        interfaceIntCoeffs,
        interfaces,
        solverDict(tolerance, relTol)
    )
{}

// ************************************************************************* //
