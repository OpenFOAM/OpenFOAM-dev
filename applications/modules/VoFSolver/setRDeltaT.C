/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2022-2024 OpenFOAM Foundation
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

#include "VoFSolver.H"
#include "fvcSurfaceIntegrate.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::solvers::VoFSolver::setRDeltaT()
{
    volScalarField& rDeltaT = trDeltaT.ref();

    const dictionary& pimpleDict = pimple.dict();

    const scalar maxCo
    (
        pimpleDict.lookupOrDefault<scalar>("maxCo", 0.9)
    );

    const volScalarField rDeltaT0("rDeltaT0", rDeltaT);

    // Set the reciprocal time-step from the local Courant number
    rDeltaT.internalFieldRef() =
        fvc::surfaceSum(mag(phi))()()/((2*maxCo)*mesh.V());

    // Clip to user-defined maximum and minimum time-steps
    scalar minRDeltaT = gMin(rDeltaT.primitiveField());
    if (pimpleDict.found("maxDeltaT") || minRDeltaT < rootVSmall)
    {
        const scalar clipRDeltaT = 1/pimpleDict.lookup<scalar>("maxDeltaT");
        rDeltaT.max(clipRDeltaT);
        minRDeltaT = max(minRDeltaT, clipRDeltaT);
    }
    if (pimpleDict.found("minDeltaT"))
    {
        const scalar clipRDeltaT = 1/pimpleDict.lookup<scalar>("minDeltaT");
        rDeltaT.min(clipRDeltaT);
        minRDeltaT = min(minRDeltaT, clipRDeltaT);
    }

    Info<< "Flow time scale min/max = "
        << gMin(1/rDeltaT.primitiveField()) << ", " << 1/minRDeltaT << endl;

    setInterfaceRDeltaT(rDeltaT);

    // Limit rate of change of time scale
    // - reduce as much as required
    // - only increase at a fraction of old time scale
    if
    (
        pimpleDict.found("rDeltaTDampingCoeff")
     && runTime.timeIndex() > runTime.startTimeIndex() + 1
    )
    {
        // Damping coefficient (1-0)
        const scalar rDeltaTDampingCoeff
        (
            pimpleDict.lookup<scalar>("rDeltaTDampingCoeff")
        );

        rDeltaT = max
        (
            rDeltaT,
            (scalar(1) - rDeltaTDampingCoeff)*rDeltaT0
        );

        Info<< "Damped flow time scale min/max = "
            << gMin(1/rDeltaT.primitiveField())
            << ", " << gMax(1/rDeltaT.primitiveField()) << endl;
    }
}


// ************************************************************************* //
