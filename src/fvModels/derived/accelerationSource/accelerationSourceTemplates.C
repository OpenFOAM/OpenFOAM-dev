/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2018-2021 OpenFOAM Foundation
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

template<class AlphaRhoFieldType>
void Foam::fv::accelerationSource::add
(
    const AlphaRhoFieldType& alphaRho,
    fvMatrix<vector>& eqn,
    const word& fieldName
) const
{
    const DimensionedField<scalar, volMesh>& V = mesh().V();

    const scalar t = mesh().time().value();
    const scalar dt = mesh().time().deltaTValue();
    const vector dU = velocity_->value(t) - velocity_->value(t - dt);
    const vector a = dU/mesh().time().deltaTValue();

    const labelList& cells = set_.cells();

    forAll(cells, i)
    {
        const label celli = cells[i];
        eqn.source()[celli] -= V[celli]*alphaRho[celli]*a;
    }
}


// ************************************************************************* //
