/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2012-2023 OpenFOAM Foundation
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

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class RhoFieldType>
void Foam::porosityModels::powerLaw::apply
(
    scalarField& Udiag,
    const scalarField& V,
    const RhoFieldType& rho,
    const vectorField& U
) const
{
    const labelList& cells = mesh_.cellZones()[zoneName_];

    const scalar C0 = C0_;
    const scalar C1m1b2 = (C1_ - 1.0)/2.0;

    forAll(cells, i)
    {
        const label celli = cells[i];

        Udiag[celli] +=
            V[celli]*rho[celli]*C0*pow(magSqr(U[celli]), C1m1b2);
    }
}


template<class RhoFieldType>
void Foam::porosityModels::powerLaw::apply
(
    tensorField& AU,
    const RhoFieldType& rho,
    const vectorField& U
) const
{
    const labelList& cells = mesh_.cellZones()[zoneName_];

    const scalar C0 = C0_;
    const scalar C1m1b2 = (C1_ - 1.0)/2.0;

    forAll(cells, i)
    {
        const label celli = cells[i];

        AU[celli] =
            AU[celli] + I*(rho[celli]*C0*pow(magSqr(U[celli]), C1m1b2));
    }
}


// ************************************************************************* //
