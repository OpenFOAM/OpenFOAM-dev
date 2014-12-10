/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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

#include "anisotropicFilter.H"
#include "addToRunTimeSelectionTable.H"
#include "zeroGradientFvPatchFields.H"
#include "wallFvPatch.H"
#include "fvc.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(anisotropicFilter, 0);
    addToRunTimeSelectionTable(LESfilter, anisotropicFilter, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::anisotropicFilter::anisotropicFilter
(
    const fvMesh& mesh,
    scalar widthCoeff
)
:
    LESfilter(mesh),
    widthCoeff_(widthCoeff),
    coeff_
    (
        IOobject
        (
            "anisotropicFilterCoeff",
            mesh.time().timeName(),
            mesh
        ),
        mesh,
        dimensionedVector("zero", dimLength*dimLength, vector::zero),
        calculatedFvPatchVectorField::typeName
    )
{
    for (direction d=0; d<vector::nComponents; d++)
    {
        coeff_.internalField().replace
        (
            d,
            (1/widthCoeff_)*
            sqr
            (
                2.0*mesh.V()
               /fvc::surfaceSum(mag(mesh.Sf().component(d)))().internalField()
            )
        );
    }
}


Foam::anisotropicFilter::anisotropicFilter
(
    const fvMesh& mesh,
    const dictionary& bd
)
:
    LESfilter(mesh),
    widthCoeff_(readScalar(bd.subDict(type() + "Coeffs").lookup("widthCoeff"))),
    coeff_
    (
        IOobject
        (
            "anisotropicFilterCoeff",
            mesh.time().timeName(),
            mesh
        ),
        mesh,
        dimensionedVector("zero", dimLength*dimLength, vector::zero),
        calculatedFvPatchScalarField::typeName
    )
{
    for (direction d=0; d<vector::nComponents; d++)
    {
        coeff_.internalField().replace
        (
            d,
            (1/widthCoeff_)*
            sqr
            (
                2.0*mesh.V()
               /fvc::surfaceSum(mag(mesh.Sf().component(d)))().internalField()
            )
        );
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::anisotropicFilter::read(const dictionary& bd)
{
    bd.subDict(type() + "Coeffs").lookup("widthCoeff") >> widthCoeff_;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::anisotropicFilter::operator()
(
    const tmp<volScalarField>& unFilteredField
) const
{
    tmp<volScalarField> tmpFilteredField =
        unFilteredField
      + (
           coeff_
         & fvc::surfaceIntegrate
           (
               mesh().Sf()
              *fvc::snGrad(unFilteredField())
           )
        );

    unFilteredField.clear();

    return tmpFilteredField;
}


Foam::tmp<Foam::volVectorField> Foam::anisotropicFilter::operator()
(
    const tmp<volVectorField>& unFilteredField
) const
{
    tmp<volVectorField> tmpFilteredField =
        unFilteredField
      + (
           coeff_
         & fvc::surfaceIntegrate
           (
               mesh().Sf()
              *fvc::snGrad(unFilteredField())
           )
        );

    unFilteredField.clear();

    return tmpFilteredField;
}


Foam::tmp<Foam::volSymmTensorField> Foam::anisotropicFilter::operator()
(
    const tmp<volSymmTensorField>& unFilteredField
) const
{
    tmp<volSymmTensorField> tmpFilteredField
    (
        new volSymmTensorField
        (
            IOobject
            (
                "anisotropicFilteredSymmTensorField",
                mesh().time().timeName(),
                mesh()
            ),
            mesh(),
            unFilteredField().dimensions()
        )
    );

    for (direction d=0; d<symmTensor::nComponents; d++)
    {
        tmpFilteredField().replace
        (
            d, anisotropicFilter::operator()(unFilteredField().component(d))
        );
    }

    unFilteredField.clear();

    return tmpFilteredField;
}


Foam::tmp<Foam::volTensorField> Foam::anisotropicFilter::operator()
(
    const tmp<volTensorField>& unFilteredField
) const
{
    tmp<volTensorField> tmpFilteredField
    (
        new volTensorField
        (
            IOobject
            (
                "anisotropicFilteredTensorField",
                mesh().time().timeName(),
                mesh()
            ),
            mesh(),
            unFilteredField().dimensions()
        )
    );

    for (direction d=0; d<tensor::nComponents; d++)
    {
        tmpFilteredField().replace
        (
            d, anisotropicFilter::operator()(unFilteredField().component(d))
        );
    }

    unFilteredField.clear();

    return tmpFilteredField;
}


// ************************************************************************* //
