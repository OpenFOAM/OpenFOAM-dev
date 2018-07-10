/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2013-2018 OpenFOAM Foundation
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

#include "turbulentBreakUp.H"
#include "fvmSup.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace diameterModels
{
namespace IATEsources
{
    defineTypeNameAndDebug(turbulentBreakUp, 0);
    addToRunTimeSelectionTable(IATEsource, turbulentBreakUp, dictionary);
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::diameterModels::IATEsources::turbulentBreakUp::
turbulentBreakUp
(
    const IATE& iate,
    const dictionary& dict
)
:
    IATEsource(iate),
    Cti_("Cti", dimless, dict),
    WeCr_("WeCr", dimless, dict)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::fvScalarMatrix>
Foam::diameterModels::IATEsources::turbulentBreakUp::R
(
    const volScalarField& alphai,
    volScalarField& kappai
) const
{
    volScalarField::Internal R
    (
        IOobject
        (
            "turbulentBreakUp:R",
            iate_.phase().time().timeName(),
            iate_.phase().mesh()
        ),
        iate_.phase().mesh(),
        dimensionedScalar("R", kappai.dimensions()/dimTime, 0)
    );

    const scalar Cti = Cti_.value();
    const scalar WeCr = WeCr_.value();
    const volScalarField Ut(this->Ut());
    const volScalarField We(this->We());

    forAll(R, celli)
    {
        if (We[celli] > WeCr)
        {
            R[celli] =
               2*Cti*Ut[celli]*sqrt(1 - WeCr/We[celli])*exp(-WeCr/We[celli]);
        }
    }

    return fvm::Su(R, kappai);
}


// ************************************************************************* //
