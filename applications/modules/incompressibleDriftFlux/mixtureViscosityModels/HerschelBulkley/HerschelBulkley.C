/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2024 OpenFOAM Foundation
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

#include "HerschelBulkley.H"
#include "incompressibleDriftFluxMixture.H"
#include "fvcGrad.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace mixtureViscosityModels
{
    defineTypeNameAndDebug(HerschelBulkley, 0);

    addToRunTimeSelectionTable
    (
        mixtureViscosityModel,
        HerschelBulkley,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::mixtureViscosityModels::HerschelBulkley::HerschelBulkley
(
    const incompressibleDriftFluxMixture& mixture
)
:
    mixtureViscosityModel(mixture),
    n_("n", dimless, coeffDict()),
    k_("k", dimDynamicViscosity*pow(dimTime, n_ - 1), coeffDict()),
    tau0_("tau0", dimDynamicViscosity/dimTime, coeffDict())
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::mixtureViscosityModels::HerschelBulkley::mu
(
    const volScalarField& muc,
    const volVectorField& U
) const
{
    const volScalarField strainRate(sqrt(2.0)*mag(symm(fvc::grad(U))));

    return min
    (
        muc,
        (tau0_ + k_*pow(strainRate, n_))
       /max(strainRate, dimensionedScalar(dimless/dimTime, rootVSmall))
    );
}


bool Foam::mixtureViscosityModels::HerschelBulkley::read()
{
    if (mixtureViscosityModel::read())
    {
        const dictionary& dict = coeffDict();

        n_.read(dict);
        k_.read(dict);
        tau0_.read(dict);

        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
