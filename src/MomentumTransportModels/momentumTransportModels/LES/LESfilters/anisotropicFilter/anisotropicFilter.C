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

#include "anisotropicFilter.H"
#include "fviSurfaceIntegrate.H"
#include "fvcSurfaceIntegrate.H"
#include "fvcSnGrad.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(anisotropicFilter, 0);
    addToRunTimeSelectionTable(LESfilter, anisotropicFilter, dictionary);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
inline Foam::tmp<Foam::VolInternalField<Type>> Foam::anisotropicFilter::filter
(
    const VolField<Type>& unFilteredField
) const
{
    tmp<VolInternalField<Type>> tmpFilteredField
    (
        VolInternalField<Type>::New
        (
            "anisotropicFilteredField",
            mesh(),
            unFilteredField.dimensions()
        )
    );

    for (direction d=0; d<pTraits<Type>::nComponents; d++)
    {
        tmpFilteredField.ref().replace
        (
            d, anisotropicFilter::operator[](unFilteredField.component(d))
        );
    }

    return tmpFilteredField;
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
            mesh.time().name(),
            mesh
        ),
        mesh,
        dimensionedVector(dimensions::area, Zero)
    )
{
    for (direction d=0; d<vector::nComponents; d++)
    {
        coeff_.primitiveFieldRef().replace
        (
            d,
            (1/widthCoeff_)*
            sqr
            (
                2*mesh.V().primitiveField()
               /fvi::surfaceSum(mag(mesh.Sf().component(d)))().primitiveField()
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
    widthCoeff_
    (
        bd.optionalTypeDict(type()).lookup<scalar>("widthCoeff")
    ),
    coeff_
    (
        IOobject
        (
            "anisotropicFilterCoeff",
            mesh.time().name(),
            mesh
        ),
        mesh,
        dimensionedVector(dimensions::area, Zero)
    )
{
    for (direction d=0; d<vector::nComponents; d++)
    {
        coeff_.replace
        (
            d,
            (1/widthCoeff_)
           *sqr
            (
                2*mesh.V()
               /fvi::surfaceSum(mag(mesh.Sf().component(d)))
            )
        );
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::anisotropicFilter::read(const dictionary& bd)
{
    bd.optionalTypeDict(type()).lookup("widthCoeff") >> widthCoeff_;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

Foam::tmp<Foam::volInternalScalarField> Foam::anisotropicFilter::operator[]
(
    const volScalarField& unFilteredField
) const
{
    return
        unFilteredField()
      + (
           coeff_
         & fvi::surfaceIntegrate
           (
               mesh().Sf()*fvc::snGrad(unFilteredField)
           )
        );
}


Foam::tmp<Foam::volInternalVectorField> Foam::anisotropicFilter::operator[]
(
    const volVectorField& unFilteredField
) const
{
    return filter(unFilteredField);
}


Foam::tmp<Foam::volInternalSymmTensorField> Foam::anisotropicFilter::operator[]
(
    const volSymmTensorField& unFilteredField
) const
{
    return filter(unFilteredField);
}


Foam::tmp<Foam::volInternalTensorField> Foam::anisotropicFilter::operator[]
(
    const volTensorField& unFilteredField
) const
{
    return filter(unFilteredField);
}


// ************************************************************************* //
