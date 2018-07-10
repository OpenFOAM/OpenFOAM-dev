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

#include "exponential.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace ParticleStressModels
{
    defineTypeNameAndDebug(exponential, 0);

    addToRunTimeSelectionTable
    (
        ParticleStressModel,
        exponential,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::ParticleStressModels::exponential::exponential
(
    const dictionary& dict
)
:
    ParticleStressModel(dict),
    preExp_(readScalar(dict.lookup("preExp"))),
    expMax_(readScalar(dict.lookup("expMax"))),
    g0_(readScalar(dict.lookup("g0")))
{}


Foam::ParticleStressModels::exponential::exponential
(
    const exponential& hc
)
:
    ParticleStressModel(hc),
    preExp_(hc.preExp_),
    expMax_(hc.expMax_),
    g0_(hc.g0_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::ParticleStressModels::exponential::~exponential()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::Field<Foam::scalar>>
Foam::ParticleStressModels::exponential::tau
(
    const Field<scalar>& alpha,
    const Field<scalar>& rho,
    const Field<scalar>& uSqr
) const
{
    return dTaudTheta(alpha, rho, uSqr)/preExp_;
}


Foam::tmp<Foam::Field<Foam::scalar>>
Foam::ParticleStressModels::exponential::dTaudTheta
(
    const Field<scalar>& alpha,
    const Field<scalar>& rho,
    const Field<scalar>& uSqr
) const
{
    return
        g0_
       *min
        (
            exp(preExp_*(alpha - alphaPacked_)),
            expMax_
        );
}


// ************************************************************************* //
