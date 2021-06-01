/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2021 OpenFOAM Foundation
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

Class
    Foam::hydrostaticInitialisation

Description
    Optional hydrostatic initialisation of p_rgh and p by solving for and
    caching the hydrostatic ph_rgh and updating the density such that

        p = ph_rgh + rho*gh + pRef

    This initialisation process is applied at the beginning of the run (not on
    restart) if the \c hydrostaticInitialisation switch is set true in
    fvSolution/PIMPLE or fvSolution/SIMPLE.  The calculation is iterative if the
    density is a function of pressure and an optional number of iterations \c
    nHydrostaticCorrectors may be specified which defaults to 5.

SourceFiles
    hydrostaticInitialisation.C

\*---------------------------------------------------------------------------*/

#ifndef hydrostaticInitialisation_H
#define hydrostaticInitialisation_H

#include "volFieldsFwd.H"
#include "surfaceFieldsFwd.H"
#include "uniformDimensionedFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    class fluidThermo;

    void hydrostaticInitialisation
    (
        volScalarField& p_rgh,
        volScalarField& p,
        volScalarField& rho,
        const volVectorField& U,
        const volScalarField& gh,
        const surfaceScalarField& ghf,
        const uniformDimensionedScalarField& pRef,
        fluidThermo& thermo,
        const dictionary& dict
    );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
