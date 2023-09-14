/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2013-2023 OpenFOAM Foundation
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

#include "IATErandomCoalescence.H"
#include "fvmSup.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace diameterModels
{
namespace IATEsources
{
    defineTypeNameAndDebug(randomCoalescence, 0);
    addToRunTimeSelectionTable(IATEsource, randomCoalescence, dictionary);
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::diameterModels::IATEsources::randomCoalescence::
randomCoalescence
(
    const IATE& iate,
    const dictionary& dict
)
:
    IATEsource(iate),
    Crc_("Crc", dimless, dict),
    C_("C", dimless, dict),
    alphaMax_("alphaMax", dimless, dict)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::fvScalarMatrix>
Foam::diameterModels::IATEsources::randomCoalescence::R
(
    const volScalarField& alphai,
    volScalarField& kappai
) const
{
    volScalarField::Internal R
    (
        IOobject
        (
            typedName("R"),
            iate_.phase().time().name(),
            iate_.phase().mesh()
        ),
        iate_.phase().mesh(),
        dimensionedScalar(dimless/dimTime, 0)
    );

    const scalar Crc = Crc_.value();
    const scalar C = C_.value();
    const scalar alphaMax = alphaMax_.value();
    const volScalarField Ut(this->Ut());
    const volScalarField& alpha = phase();
    const scalar cbrtAlphaMax = cbrt(alphaMax);

    forAll(R, celli)
    {
        if (kappai[celli] > 0 && alpha[celli] < alphaMax - small)
        {
            const scalar cbrtAlphaMaxMAlpha = cbrtAlphaMax - cbrt(alpha[celli]);

            R[celli] =
                12*phi()*kappai[celli]*alpha[celli]
               *Crc
               *Ut[celli]
               *(1 - exp(-C*cbrt(alpha[celli]*alphaMax)/cbrtAlphaMaxMAlpha))
               /(cbrtAlphaMax*cbrtAlphaMaxMAlpha);
        }
    }

    return -fvm::Sp(R, kappai);
}


// ************************************************************************* //
