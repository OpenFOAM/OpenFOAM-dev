/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012-2013 OpenFOAM Foundation
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
void Foam::porosityModels::DarcyForchheimer::apply
(
    scalarField& Udiag,
    vectorField& Usource,
    const scalarField& V,
    const RhoFieldType& rho,
    const scalarField& mu,
    const vectorField& U
) const
{
    forAll(cellZoneIDs_, zoneI)
    {
        const tensorField& dZones = D_[zoneI];
        const tensorField& fZones = F_[zoneI];

        const labelList& cells = mesh_.cellZones()[cellZoneIDs_[zoneI]];

        forAll(cells, i)
        {
            const label cellI = cells[i];
            const label j = this->fieldIndex(i);
            const tensor Cd =
                mu[cellI]*dZones[j] + (rho[cellI]*mag(U[cellI]))*fZones[j];

            const scalar isoCd = tr(Cd);

            Udiag[cellI] += V[cellI]*isoCd;
            Usource[cellI] -= V[cellI]*((Cd - I*isoCd) & U[cellI]);
        }
    }
}


template<class RhoFieldType>
void Foam::porosityModels::DarcyForchheimer::apply
(
    tensorField& AU,
    const RhoFieldType& rho,
    const scalarField& mu,
    const vectorField& U
) const
{
    forAll(cellZoneIDs_, zoneI)
    {
        const tensorField& dZones = D_[zoneI];
        const tensorField& fZones = F_[zoneI];

        const labelList& cells = mesh_.cellZones()[cellZoneIDs_[zoneI]];

        forAll(cells, i)
        {
            const label cellI = cells[i];
            const label j = this->fieldIndex(i);
            const tensor D = dZones[j];
            const tensor F = fZones[j];

            AU[cellI] += mu[cellI]*D + (rho[cellI]*mag(U[cellI]))*F;
        }
    }
}


// ************************************************************************* //
