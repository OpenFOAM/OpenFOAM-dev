/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2026 OpenFOAM Foundation
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

#include "Gulder.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace XiEqModels
{
    defineTypeNameAndDebug(Gulder, 0);
    addToRunTimeSelectionTable(XiEqModel, Gulder, dictionary);
}
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

bool Foam::XiEqModels::Gulder::readCoeffs(const dictionary& dict)
{
    XiEqModel::readCoeffs(dict);

    XiEqCoeff_.readIfPresent(dict);
    SuMin_.readIfPresent(dict);

    return true;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::XiEqModels::Gulder::Gulder
(
    const dictionary& dict,
    const ubRhoThermo& thermo,
    const compressibleMomentumTransportModel& turbulence,
    const volScalarField& Su
)
:
    XiEqModel(thermo, turbulence, Su),
    XiEqCoeff_("XiEqCoeff", dimless, 0.62),
    SuMin_("SuMin", 0.01*Su.average())
{
    readCoeffs(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::XiEqModels::Gulder::~Gulder()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::XiEqModels::Gulder::XiEq() const
{
    const volScalarField up(sqrt((2.0/3.0)*momentumTransport_.k()));

    const volScalarField epsilon
    (
        max
        (
            momentumTransport_.epsilon(),
            dimensionedScalar(sqr(dimVelocity)/dimTime, small)
        )
    );

    const volScalarField tauEta
    (
        sqrt(thermo_.uThermo().mu()/(thermo_.uThermo().rho()*epsilon))
    );

    const volScalarField Reta(up/sqrt(epsilon*tauEta));

    return (1 + XiEqCoeff_*sqrt(up/max(Su_, SuMin_))*Reta);
}


// ************************************************************************* //
