/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2018-2020 OpenFOAM Foundation
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

scalarList W(const fluidReactionThermo& thermo)
{
    const basicSpecieMixture& composition = thermo.composition();

    scalarList W(composition.Y().size());

    forAll(W, i)
    {
        W[i] = composition.Wi(i);
    }

    return W;
}


scalar h0
(
    const fluidReactionThermo& thermo,
    const scalarList& Y,
    const scalar p,
    const scalar T
)
{
    const basicSpecieMixture& composition = thermo.composition();

    scalar h0 = 0;

    forAll(Y, i)
    {
        h0 += Y[i]*composition.Hs(i, p, T);
    }

    return h0;
}


// ************************************************************************* //
