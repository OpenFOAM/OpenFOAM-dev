/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2022 OpenFOAM Foundation
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

#include "simpleFilter.H"
#include "fvcSurfaceIntegrate.H"
#include "surfaceInterpolate.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(simpleFilter, 0);
    addToRunTimeSelectionTable(LESfilter, simpleFilter, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::simpleFilter::simpleFilter
(
    const fvMesh& mesh
)
:
    LESfilter(mesh)
{}


Foam::simpleFilter::simpleFilter(const fvMesh& mesh, const dictionary&)
:
    LESfilter(mesh)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::simpleFilter::read(const dictionary&)
{}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::simpleFilter::operator()
(
    const tmp<volScalarField>& unFilteredField
) const
{
    correctBoundaryConditions(unFilteredField);

    tmp<volScalarField> filteredField = fvc::surfaceSum
    (
        mesh().magSf()*fvc::interpolate(unFilteredField)
    )/fvc::surfaceSum(mesh().magSf());

    unFilteredField.clear();

    return filteredField;
}


Foam::tmp<Foam::volVectorField> Foam::simpleFilter::operator()
(
    const tmp<volVectorField>& unFilteredField
) const
{
    correctBoundaryConditions(unFilteredField);

    tmp<volVectorField> filteredField = fvc::surfaceSum
    (
        mesh().magSf()*fvc::interpolate(unFilteredField)
    )/fvc::surfaceSum(mesh().magSf());

    unFilteredField.clear();

    return filteredField;
}


Foam::tmp<Foam::volSymmTensorField> Foam::simpleFilter::operator()
(
    const tmp<volSymmTensorField>& unFilteredField
) const
{
    correctBoundaryConditions(unFilteredField);

    tmp<volSymmTensorField> filteredField = fvc::surfaceSum
    (
        mesh().magSf()*fvc::interpolate(unFilteredField)
    )/fvc::surfaceSum(mesh().magSf());

    unFilteredField.clear();

    return filteredField;
}


Foam::tmp<Foam::volTensorField> Foam::simpleFilter::operator()
(
    const tmp<volTensorField>& unFilteredField
) const
{
    correctBoundaryConditions(unFilteredField);

    tmp<volTensorField> filteredField = fvc::surfaceSum
    (
        mesh().magSf()*fvc::interpolate(unFilteredField)
    )/fvc::surfaceSum(mesh().magSf());

    unFilteredField.clear();

    return filteredField;
}


// ************************************************************************* //
