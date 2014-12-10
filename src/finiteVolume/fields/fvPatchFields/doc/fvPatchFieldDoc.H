/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012 OpenFOAM Foundation
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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

\defgroup grpBoundaryConditions Boundary Conditions
@{
    This group contains OpenFOAM boundary condition types.  All conditions
    are derived from the base Foam::fvPatchField class.  Patch values are
    typically updated using the virtual functions
    - Foam::fvPatchField::updateCoeffs()
    - Foam::fvPatchField::evaluate()

    The difference between the methods is based on when the patch values are
    updated.  When the condition is applied to a solution variable, the
    call to updateCoeffs() occurs as a preliminary step of the
    \<matrix\>.solve().  The evaluate() method is invoked after, or independent
    of the matrix solve, via a call to \<field\>.correctBoundaryConditions().
@}

\defgroup grpConstraintBoundaryConditions Constraint boundary Conditions
@{
    \ingroup grpBoundaryConditions
    This group contains constraint boundary condition types.  These conditions
    are mainly employed to reduced dimensioned cases.
@}

\defgroup grpInletBoundaryConditions Inlet boundary Conditions
@{
    \ingroup grpBoundaryConditions
    This group contains inlet boundary condition types
@}

\defgroup grpOutletBoundaryConditions Outlet boundary Conditions
@{
    \ingroup grpBoundaryConditions
    This group contains outlet boundary condition types
@}

\defgroup grpGenericBoundaryConditions Generic boundary Conditions
@{
    \ingroup grpBoundaryConditions
    This group contains generic boundary condition types
@}

\defgroup grpCoupledBoundaryConditions Coupled boundary Conditions
@{
    \ingroup grpBoundaryConditions
    This group contains coupled boundary condition types
@}

\defgroup grpWallBoundaryConditions Wall boundary Conditions
@{
    \ingroup grpBoundaryConditions
    This group contains wall boundary condition types
@}

\*---------------------------------------------------------------------------*/
