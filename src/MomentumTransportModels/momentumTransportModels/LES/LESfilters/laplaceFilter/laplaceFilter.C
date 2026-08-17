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

\*---------------------------------------------------------------------------*/

#include "laplaceFilter.H"
#include "fviLaplacian.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(laplaceFilter, 0);
    addToRunTimeSelectionTable(LESfilter, laplaceFilter, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::laplaceFilter::laplaceFilter(const fvMesh& mesh, scalar widthCoeff)
:
    LESfilter(mesh),
    widthCoeff_(widthCoeff),
    coeff_
    (
        IOobject
        (
            "laplaceFilterCoeff",
            mesh.time().name(),
            mesh
        ),
        mesh,
        dimensionedScalar(dimensions::area, 0),
        calculatedFvPatchScalarField::typeName
    )
{
    coeff_.internalFieldRef() = pow(mesh.V(), 2.0/3.0)/widthCoeff_;
}


Foam::laplaceFilter::laplaceFilter(const fvMesh& mesh, const dictionary& bd)
:
    LESfilter(mesh),
    widthCoeff_
    (
        bd.optionalTypeDict(type()).lookup<scalar>("widthCoeff")
    ),
    coeff_
    (
        IOobject
        (
            "laplaceFilterCoeff",
            mesh.time().name(),
            mesh
        ),
        mesh,
        dimensionedScalar(dimensions::area, 0),
        calculatedFvPatchScalarField::typeName
    )
{
    coeff_.internalFieldRef() = pow(mesh.V(), 2.0/3.0)/widthCoeff_;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::laplaceFilter::read(const dictionary& bd)
{
    bd.optionalTypeDict(type()).lookup("widthCoeff") >> widthCoeff_;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

Foam::tmp<Foam::volInternalScalarField> Foam::laplaceFilter::operator[]
(
    const volScalarField& unFilteredField
) const
{
    return unFilteredField() + fvi::laplacian(coeff_, unFilteredField);
}


Foam::tmp<Foam::volInternalVectorField> Foam::laplaceFilter::operator[]
(
    const volVectorField& unFilteredField
) const
{
    return unFilteredField() + fvi::laplacian(coeff_, unFilteredField);
}


Foam::tmp<Foam::volInternalSymmTensorField> Foam::laplaceFilter::operator[]
(
    const volSymmTensorField& unFilteredField
) const
{
    return unFilteredField() + fvi::laplacian(coeff_, unFilteredField);
}


Foam::tmp<Foam::volInternalTensorField> Foam::laplaceFilter::operator[]
(
    const volTensorField& unFilteredField
) const
{
    return unFilteredField() + fvi::laplacian(coeff_, unFilteredField);
}


// ************************************************************************* //
