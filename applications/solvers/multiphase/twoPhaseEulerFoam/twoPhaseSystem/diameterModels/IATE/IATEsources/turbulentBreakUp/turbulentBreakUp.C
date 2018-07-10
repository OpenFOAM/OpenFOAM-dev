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

Foam::tmp<Foam::volScalarField>
Foam::diameterModels::IATEsources::turbulentBreakUp::R() const
{
    tmp<volScalarField> tR
    (
        new volScalarField
        (
            IOobject
            (
                "R",
                iate_.phase().U().time().timeName(),
                iate_.phase().mesh()
            ),
            iate_.phase().U().mesh(),
            dimensionedScalar("R", dimless/dimTime, 0)
        )
    );

    volScalarField R = tR();

    scalar Cti = Cti_.value();
    scalar WeCr = WeCr_.value();
    volScalarField Ut(this->Ut());
    volScalarField We(this->We());
    const volScalarField& d(iate_.d()());

    forAll(R, celli)
    {
        if (We[celli] > WeCr)
        {
            R[celli] =
                (1.0/3.0)
               *Cti/d[celli]
               *Ut[celli]
               *sqrt(1 - WeCr/We[celli])
               *exp(-WeCr/We[celli]);
        }
    }

    return tR;
}


// ************************************************************************* //
