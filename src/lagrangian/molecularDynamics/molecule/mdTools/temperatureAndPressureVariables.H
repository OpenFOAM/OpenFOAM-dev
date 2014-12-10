/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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
    temperatureAndPressureVariables.H

Description
    Provides accumulation variables for temperatureAndPressure.H

\*---------------------------------------------------------------------------*/

vector accumulatedTotalLinearMomentum(vector::zero);

scalar accumulatedTotalMass = 0.0;

scalar accumulatedTotalAngularKE = 0.0;

scalar accumulatedTotalLinearKE = 0.0;

scalar accumulatedTotalPE = 0.0;

scalar accumulatedTotalrDotfSum = 0.0;

label accumulatedNMols = 0;

label accumulatedDOFs = 0;

scalar averageTemperature = 0.0;

scalar averagePressure = 0.0;

const scalarField& cellVols = mesh.cellVolumes();

scalar meshVolume = sum(cellVols);

if (Pstream::parRun())
{
    reduce(meshVolume, sumOp<scalar>());
}


// ************************************************************************* //
