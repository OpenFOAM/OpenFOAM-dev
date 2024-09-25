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

#include "Quemada.H"
#include "incompressibleDriftFluxMixture.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace mixtureViscosityModels
{
    defineTypeNameAndDebug(Quemada, 0);

    addToRunTimeSelectionTable
    (
        mixtureViscosityModel,
        Quemada,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::mixtureViscosityModels::Quemada::Quemada
(
    const incompressibleDriftFluxMixture& mixture
)
:
    mixtureViscosityModel(mixture),
    q_(optionalSubDict(typeName + "Coeffs").lookupOrDefault("q", scalar(2))),
    muMax_
    (
        "muMax",
        dimDynamicViscosity,
        optionalSubDict(typeName + "Coeffs").lookup("muMax")
    )
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::mixtureViscosityModels::Quemada::mu
(
    const volScalarField& muc,
    const volVectorField& U
) const
{
    return min
    (
        muc*pow(max(1 - mixture_.alphad()/mixture_.alphaMax(), small), -q_),
        muMax_
    );
}


bool Foam::mixtureViscosityModels::Quemada::read()
{
    if (mixtureViscosityModel::read())
    {
        const dictionary& dict = coeffDict();

        dict.lookup("q") >> q_;
        muMax_.read(dict);

        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
