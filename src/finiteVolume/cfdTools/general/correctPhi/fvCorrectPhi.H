/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2015-2023 OpenFOAM Foundation
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

Function
    Foam::fv::correctPhi

Description
    Flux correction functions to ensure continuity.

    Required during start-up, restart, mesh-motion etc. when non-conservative
    fluxes may adversely affect the prediction-part of the solution algorithm
    (the part before the first pressure solution which would ensure continuity).
    This is particularly important for VoF and other multi-phase solver in
    which non-conservative fluxes cause unboundedness of the phase-fraction.

SourceFiles
    fvCorrectPhi.C
    correctUphiBCs.C

\*---------------------------------------------------------------------------*/

#ifndef fvCorrectPhi_H
#define fvCorrectPhi_H

#include "volFieldsFwd.H"
#include "surfaceFieldsFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    class pressureReference;
    class nonOrthogonalSolutionControl;

    //- If the mesh is moving correct the velocity BCs on the moving walls to
    //  ensure the corrected fluxes and velocity are consistent
    void correctUphiBCs
    (
        volVectorField& U,
        surfaceScalarField& phi,
        const bool evaluateUBCs
    );

    //- If the mesh is moving correct the velocity BCs on the moving walls to
    //  ensure the corrected fluxes and velocity are consistent
    void correctUphiBCs
    (
        const volScalarField& rho,
        volVectorField& U,
        surfaceScalarField& phi,
        const bool evaluateUBCs
    );

    namespace fv
    {
        void correctPhi
        (
            surfaceScalarField& phi,
            const volVectorField& U,
            const volScalarField& p,
            const autoPtr<volScalarField>& rAU,
            const autoPtr<volScalarField>& divU,
            const pressureReference& pressureReference,
            nonOrthogonalSolutionControl& pcorrControl
        );

        void correctPhi
        (
            surfaceScalarField& phi,
            const volScalarField& p,
            const volScalarField& psi,
            const autoPtr<volScalarField>& rAU,
            const volScalarField& divRhoU,
            nonOrthogonalSolutionControl& pcorrControl
        );
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
