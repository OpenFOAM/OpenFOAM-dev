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

#include "TomiyamaLift.H"
#include "aspectRatioModel.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace liftModels
{
    defineTypeNameAndDebug(TomiyamaLift, 0);
    addToRunTimeSelectionTable(liftModel, TomiyamaLift, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::liftModels::TomiyamaLift::TomiyamaLift
(
    const dictionary& dict,
    const phaseInterface& interface
)
:
    dispersedLiftModel(dict, interface),
    aspectRatio_(aspectRatioModel::New(dict.subDict("aspectRatio"), interface))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::liftModels::TomiyamaLift::~TomiyamaLift()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::liftModels::TomiyamaLift::Cl() const
{
    const volScalarField EoH
    (
        interface_.Eo(interface_.dispersed().d()/cbrt(aspectRatio_->E()))
    );

    const volScalarField f
    (
        0.0010422*pow3(EoH) - 0.0159*sqr(EoH) - 0.0204*EoH + 0.474
    );

    return
        neg(EoH - 4)*min(0.288*tanh(0.121*interface_.Re()), f)
      + pos0(EoH - 4)*neg(EoH - 10.7)*f
      + pos0(EoH - 10.7)*(-0.288);
}


// ************************************************************************* //
