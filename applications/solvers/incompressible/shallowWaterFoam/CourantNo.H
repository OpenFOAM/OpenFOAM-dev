/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2022 OpenFOAM Foundation
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

Global
    CourantNo

Description
    Calculates and outputs the maximum Courant Number.

\*---------------------------------------------------------------------------*/

{
    const scalarField sumPhi
    (
        fvc::surfaceSum(mag(phi))().primitiveField()
       /h.primitiveField()
    );

    const scalar CoNum =
        0.5*gMax(sumPhi/mesh.V().field())*runTime.deltaTValue();

    const scalar meanCoNum =
        0.5*(gSum(sumPhi)/gSum(mesh.V().field()))*runTime.deltaTValue();

    // Gravity wave Courant number
    const scalar waveCoNum = 0.25*gMax
    (
        fvc::surfaceSum
        (
            fvc::interpolate(sqrt(h))*mesh.magSf()
        )().primitiveField()/mesh.V().field()
    )*sqrt(magg).value()*runTime.deltaTValue();

    Info<< "Courant number mean: " << meanCoNum
        << " max: " << CoNum << endl;

    Info<< "Gravity wave Courant number max: " << waveCoNum
        << endl;
}

// ************************************************************************* //
