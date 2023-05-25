/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2022-2023 OpenFOAM Foundation
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

#include "compressibleVoF.H"
#include "fvcDiv.H"

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::solvers::compressibleVoF::alphaSuSp
(
    tmp<volScalarField::Internal>& Su,
    tmp<volScalarField::Internal>& Sp
)
{
    Sp = volScalarField::Internal::New
    (
        "Sp",
        mesh,
        dimensionedScalar(dgdt.dimensions(), 0)
    );

    Su = volScalarField::Internal::New
    (
        "Su",
        mesh,
        dimensionedScalar(dgdt.dimensions(), 0)
    );

    if (fvModels().addsSupToField(alpha1.name()))
    {
        // Phase change alpha1 source
        const fvScalarMatrix alphaSup(fvModels().source(alpha1));

        Su = alphaSup.Su();
        Sp = alphaSup.Sp();
    }

    volScalarField::Internal& SpRef = Sp.ref();
    volScalarField::Internal& SuRef = Su.ref();

    forAll(dgdt, celli)
    {
        if (dgdt[celli] > 0.0)
        {
            SpRef[celli] -= dgdt[celli]/max(1.0 - alpha1[celli], 1e-4);
            SuRef[celli] += dgdt[celli]/max(1.0 - alpha1[celli], 1e-4);
        }
        else if (dgdt[celli] < 0.0)
        {
            SpRef[celli] += dgdt[celli]/max(alpha1[celli], 1e-4);
        }
    }
}


// ************************************************************************* //
