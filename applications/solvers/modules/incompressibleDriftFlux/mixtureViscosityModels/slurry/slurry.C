/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2014-2023 OpenFOAM Foundation
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

#include "slurry.H"
#include "incompressibleDriftFluxMixture.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace mixtureViscosityModels
{
    defineTypeNameAndDebug(slurry, 0);

    addToRunTimeSelectionTable
    (
        mixtureViscosityModel,
        slurry,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::mixtureViscosityModels::slurry::slurry
(
    const incompressibleDriftFluxMixture& mixture
)
:
    mixtureViscosityModel(mixture)
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::mixtureViscosityModels::slurry::mu
(
    const volScalarField& muc,
    const volVectorField& U
) const
{
    const volScalarField& alphad = mixture_.alphad();

    return
    (
        muc*(1.0 + 2.5*alphad + 10.05*sqr(alphad) + 0.00273*exp(16.6*alphad))
    );
}


bool Foam::mixtureViscosityModels::slurry::read()
{
    return mixtureViscosityModel::read();
}


// ************************************************************************* //
