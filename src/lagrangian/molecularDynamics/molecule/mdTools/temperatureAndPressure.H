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
    temperatureAndPressure.H

Description
    Accumulates values for temperature and pressure measurement, and
    calculates and outputs the average values at output times.
    Requires temperatureAndPressureVariables.H to be declared before the
    timeloop.

\*---------------------------------------------------------------------------*/

accumulatedTotalLinearMomentum += singleStepTotalLinearMomentum;

accumulatedTotalMass += singleStepTotalMass;

accumulatedTotalLinearKE += singleStepTotalLinearKE;

accumulatedTotalAngularKE += singleStepTotalAngularKE;

accumulatedTotalPE += singleStepTotalPE;

accumulatedTotalrDotfSum += singleStepTotalrDotf;

accumulatedNMols += singleStepNMols;

accumulatedDOFs += singleStepDOFs;

if (runTime.outputTime())
{
    if (accumulatedNMols)
    {
        Info<< "calculating averages" << endl;

        averageTemperature =
        (
            2.0/(physicoChemical::k.value()*accumulatedDOFs)
            *
            (
                accumulatedTotalLinearKE + accumulatedTotalAngularKE
                -
                0.5*magSqr(accumulatedTotalLinearMomentum)/accumulatedTotalMass
            )
        );

        averagePressure =
        (
            (
                (accumulatedNMols/nAveragingSteps)
               *physicoChemical::k.value()*averageTemperature
              + accumulatedTotalrDotfSum/(6.0*nAveragingSteps)
            )
            /
            meshVolume
        );

        Info<< "----------------------------------------" << nl
            << "Averaged properties" << nl
            << "Average |velocity| = "
            << mag(accumulatedTotalLinearMomentum)/accumulatedTotalMass << nl
            << "Average temperature = " << averageTemperature << nl
            << "Average pressure = " << averagePressure << nl
            << "----------------------------------------" << endl;
    }
    else
    {
        Info<< "Not averaging temperature and pressure: "
            << "no molecules in system" << endl;
    }

    accumulatedTotalLinearMomentum = vector::zero;

    accumulatedTotalMass = 0.0;

    accumulatedTotalLinearKE = 0.0;

    accumulatedTotalAngularKE = 0.0;

    accumulatedTotalPE = 0.0;

    accumulatedTotalrDotfSum = 0.0;

    accumulatedNMols = 0;

    accumulatedDOFs = 0;
}


// ************************************************************************* //
