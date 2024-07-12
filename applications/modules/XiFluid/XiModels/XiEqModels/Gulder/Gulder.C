/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2024 OpenFOAM Foundation
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

    XiEqCoeff_ = dict.lookupOrDefault<scalar>("XiEqCoeff", 0.62);
    uPrimeCoeff_ = dict.lookupOrDefault<scalar>("uPrimeCoeff", 1);

    return true;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::XiEqModels::Gulder::Gulder
(
    const dictionary& dict,
    const psiuMulticomponentThermo& thermo,
    const fluidThermoThermophysicalTransportModel& turbulence,
    const volScalarField& Su
)
:
    XiEqModel(thermo, turbulence, Su),
    XiEqCoeff_(dict.lookupOrDefault<scalar>("XiEqCoeff", 0.62)),
    uPrimeCoeff_(dict.lookupOrDefault<scalar>("uPrimeCoeff", 1)),
    SuMin_(0.01*Su.average())
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::XiEqModels::Gulder::~Gulder()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::XiEqModels::Gulder::XiEq() const
{
    const volScalarField up(sqrt((2.0/3.0)*turbulence_.k()));

    tmp<volScalarField> tepsilon = turbulence_.epsilon();
    const volScalarField& epsilon = tepsilon();

    const volScalarField tauEta
    (
        sqrt(mag(thermo_.muu()/(thermo_.rhou()*epsilon)))
    );

    const volScalarField Reta
    (
        up
       /(
            sqrt(epsilon*tauEta)
          + dimensionedScalar(up.dimensions(), 1e-8)
        )
    );

    return (1 + XiEqCoeff_*sqrt(up/(Su_ + SuMin_))*Reta);
}


// ************************************************************************* //
