/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2023-2024 OpenFOAM Foundation
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

#include "multiphaseVoFSolver.H"
#include "fvcSmooth.H"
#include "fvcSurfaceIntegrate.H"
#include "fvcAverage.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::solvers::multiphaseVoFSolver::setInterfaceRDeltaT
(
    volScalarField& rDeltaT
)
{
    const dictionary& pimpleDict = pimple.dict();

    const scalar maxCo
    (
        pimpleDict.lookupOrDefault<scalar>("maxCo", 0.9)
    );

    const scalar maxAlphaCo
    (
        pimpleDict.lookupOrDefault<scalar>("maxAlphaCo", 0.2)
    );

    const scalar rDeltaTSmoothingCoeff
    (
        pimpleDict.lookupOrDefault<scalar>("rDeltaTSmoothingCoeff", 0.1)
    );

    const scalar alphaSpreadMax
    (
        pimpleDict.lookupOrDefault<scalar>("alphaSpreadMax", 0.99)
    );

    const scalar alphaSpreadMin
    (
        pimpleDict.lookupOrDefault<scalar>("alphaSpreadMin", 0.01)
    );

    if (maxAlphaCo < maxCo)
    {
        // Further limit the reciprocal time-step
        // in the vicinity of the interface

        const volScalarField::Internal alphaCoRdeltaT
        (
            fvc::surfaceSum(mag(phi))()()/((2*maxAlphaCo)*mesh.V())
        );

        forAll(phases, phasei)
        {
            const volScalarField alphaBar(fvc::average(phases[phasei]));

            rDeltaT.internalFieldRef() =
                max
                (
                    rDeltaT(),
                    pos0(alphaBar() - alphaSpreadMin)
                   *pos0(alphaSpreadMax - alphaBar())
                   *alphaCoRdeltaT
                );
        }
    }

    // Update the boundary values of the reciprocal time-step
    rDeltaT.correctBoundaryConditions();

    Info<< "Flow and interface time scale min/max = "
        << gMin(1/rDeltaT.primitiveField())
        << ", " << gMax(1/rDeltaT.primitiveField()) << endl;

    if (rDeltaTSmoothingCoeff < 1.0)
    {
        fvc::smooth(rDeltaT, rDeltaTSmoothingCoeff);
    }

    Info<< "Smoothed flow time scale min/max = "
        << gMin(1/rDeltaT.primitiveField())
        << ", " << gMax(1/rDeltaT.primitiveField()) << endl;
}


// ************************************************************************* //
