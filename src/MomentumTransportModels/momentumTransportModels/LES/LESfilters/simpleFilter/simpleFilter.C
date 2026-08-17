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

#include "simpleFilter.H"
#include "fviSurfaceIntegrate.H"
#include "surfaceInterpolate.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(simpleFilter, 0);
    addToRunTimeSelectionTable(LESfilter, simpleFilter, dictionary);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
inline Foam::tmp<Foam::VolInternalField<Type>> Foam::simpleFilter::filter
(
    const VolField<Type>& unFilteredField
) const
{
    tmp<VolInternalField<Type>> filteredField
    (
        VolInternalField<Type>::New
        (
            "simpleFilteredField",
            fvi::surfaceSum
            (
                mesh().magSf()*fvc::interpolate(unFilteredField)
            )/fvi::surfaceSum(mesh().magSf())
        )
    );

    return filteredField;
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

Foam::tmp<Foam::volInternalScalarField> Foam::simpleFilter::operator[]
(
    const volScalarField& unFilteredField
) const
{
    return filter(unFilteredField);
}


Foam::tmp<Foam::volInternalVectorField> Foam::simpleFilter::operator[]
(
    const volVectorField& unFilteredField
) const
{
    return filter(unFilteredField);
}


Foam::tmp<Foam::volInternalSymmTensorField> Foam::simpleFilter::operator[]
(
    const volSymmTensorField& unFilteredField
) const
{
    return filter(unFilteredField);
}


Foam::tmp<Foam::volInternalTensorField> Foam::simpleFilter::operator[]
(
    const volTensorField& unFilteredField
) const
{
    return filter(unFilteredField);
}


// ************************************************************************* //
