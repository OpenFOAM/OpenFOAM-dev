/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2014-2024 OpenFOAM Foundation
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

#include "BinghamPlastic.H"
#include "incompressibleDriftFluxMixture.H"
#include "fvcGrad.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace mixtureViscosityModels
{
    defineTypeNameAndDebug(BinghamPlastic, 0);

    addToRunTimeSelectionTable
    (
        mixtureViscosityModel,
        BinghamPlastic,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::mixtureViscosityModels::BinghamPlastic::BinghamPlastic
(
    const incompressibleDriftFluxMixture& mixture
)
:
    plastic(mixture),
    yieldStressCoeff_
    (
        "BinghamCoeff",
        dimensionSet(1, -1, -2, 0, 0),
        coeffDict()
    ),
    yieldStressExponent_
    (
        "BinghamExponent",
        dimless,
        coeffDict()
    ),
    yieldStressOffset_
    (
        "BinghamOffset",
        dimless,
        coeffDict()
    )
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::mixtureViscosityModels::BinghamPlastic::mu
(
    const volScalarField& muc,
    const volVectorField& U
) const
{
    volScalarField tauy
    (
        yieldStressCoeff_
       *(
            pow
            (
                scalar(10),
                min
                (
                    log10(vGreat),
                    yieldStressExponent_
                   *(max(mixture_.alphad(), scalar(0)) + yieldStressOffset_)
                )
            )
          - pow
            (
                scalar(10),
                yieldStressExponent_*yieldStressOffset_
            )
        )
    );

    volScalarField mup(plastic::mu(muc, U));

    dimensionedScalar tauySmall("tauySmall", tauy.dimensions(), small);

    return min
    (
        tauy
       /(
            sqrt(2.0)*mag(symm(fvc::grad(U)))
          + 1.0e-4*(tauy + tauySmall)/mup
        )
      + mup,
        muMax_
    );
}


bool Foam::mixtureViscosityModels::BinghamPlastic::read()
{
    if (plastic::read())
    {
        const dictionary& plasticCoeffs = coeffDict();

        yieldStressCoeff_.read(plasticCoeffs);
        yieldStressExponent_.read(plasticCoeffs);
        yieldStressOffset_.read(plasticCoeffs);

        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
