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

#include "plastic.H"
#include "incompressibleDriftFluxMixture.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace mixtureViscosityModels
{
    defineTypeNameAndDebug(plastic, 0);

    addToRunTimeSelectionTable
    (
        mixtureViscosityModel,
        plastic,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::mixtureViscosityModels::plastic::plastic
(
    const incompressibleDriftFluxMixture& mixture
)
:
    mixtureViscosityModel(mixture),
    plasticCoeffs_(optionalSubDict(typeName + "Coeffs")),
    plasticViscosityCoeff_
    (
        "coeff",
        dimDynamicViscosity,
        plasticCoeffs_.lookup("coeff")
    ),
    plasticViscosityExponent_
    (
        "exponent",
        dimless,
        plasticCoeffs_.lookup("exponent")
    ),
    muMax_("muMax", dimDynamicViscosity, plasticCoeffs_.lookup("muMax"))
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::mixtureViscosityModels::plastic::mu
(
    const volScalarField& muc,
    const volVectorField& U
) const
{
    return min
    (
        muc
      + plasticViscosityCoeff_
       *(
            pow
            (
                scalar(10),
                plasticViscosityExponent_*mixture_.alphad()
            ) - scalar(1)
        ),
        muMax_
    );
}


bool Foam::mixtureViscosityModels::plastic::read()
{
    if (mixtureViscosityModel::read())
    {
        plasticCoeffs_ = optionalSubDict(typeName + "Coeffs");

        plasticCoeffs_.lookup("k") >> plasticViscosityCoeff_;
        plasticCoeffs_.lookup("n") >> plasticViscosityExponent_;
        plasticCoeffs_.lookup("muMax") >> muMax_;

        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
