/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2017-2023 OpenFOAM Foundation
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

#include "volFields.H"
#include "geometricOneField.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class AlphaFieldType, class RhoFieldType>
void Foam::porosityModels::solidification::apply
(
    scalarField& Udiag,
    const scalarField& V,
    const AlphaFieldType& alpha,
    const RhoFieldType& rho,
    const volVectorField& U
) const
{
    const labelList& cells = mesh_.cellZones()[zoneName_];

    const volScalarField& T = mesh_.lookupObject<volScalarField>
    (
        IOobject::groupName(TName_, U.group())
    );

    forAll(cells, i)
    {
        const label celli = cells[i];
        Udiag[celli] +=
            V[celli]*alpha[celli]*rho[celli]*D_->value(T[celli]);
    }
}


template<class AlphaFieldType, class RhoFieldType>
void Foam::porosityModels::solidification::apply
(
    tensorField& AU,
    const AlphaFieldType& alpha,
    const RhoFieldType& rho,
    const volVectorField& U
) const
{
    const labelList& cells = mesh_.cellZones()[zoneName_];

    const volScalarField& T = mesh_.lookupObject<volScalarField>
    (
        IOobject::groupName(TName_, U.group())
    );

    forAll(cells, i)
    {
        const label celli = cells[i];
        AU[celli] +=
            tensor::I*alpha[celli]*rho[celli]*D_->value(T[celli]);
    }
}


template<class RhoFieldType>
void Foam::porosityModels::solidification::apply
(
    scalarField& Udiag,
    const scalarField& V,
    const RhoFieldType& rho,
    const volVectorField& U
) const
{
    if (alphaName_ == "none")
    {
        return apply(Udiag, V, geometricOneField(), rho, U);
    }
    else
    {
        const volScalarField& alpha = mesh_.lookupObject<volScalarField>
        (
            IOobject::groupName(alphaName_, U.group())
        );

        return apply(Udiag, V, alpha, rho, U);
    }
}


template<class RhoFieldType>
void Foam::porosityModels::solidification::apply
(
    tensorField& AU,
    const RhoFieldType& rho,
    const volVectorField& U
) const
{
    if (alphaName_ == "none")
    {
        return apply(AU, geometricOneField(), rho, U);
    }
    else
    {
        const volScalarField& alpha = mesh_.lookupObject<volScalarField>
        (
            IOobject::groupName(alphaName_, U.group())
        );

        return apply(AU, alpha, rho, U);
    }
}


// ************************************************************************* //
