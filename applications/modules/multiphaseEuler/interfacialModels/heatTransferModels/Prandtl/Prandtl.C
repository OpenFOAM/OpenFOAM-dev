/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2023 OpenFOAM Foundation
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

#include "Prandtl.H"
#include "dragModel.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace heatTransferModels
{
    defineTypeNameAndDebug(Prandtl, 0);
    addToRunTimeSelectionTable
    (
        heatTransferModel,
        Prandtl,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::heatTransferModels::Prandtl::Prandtl
(
    const dictionary& dict,
    const phaseInterface& interface,
    const bool registerObject
)
:
    heatTransferModel(dict, interface, registerObject),
    interfacePtr_(interface.clone()),
    interface_(interfacePtr_()),
    Pr_("Pr", dimless, dict)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::heatTransferModels::Prandtl::~Prandtl()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::heatTransferModels::Prandtl::K
(
    const scalar residualAlpha
) const
{
    const dragModel& drag =
        interface_.fluid().lookupInterfacialModel<dragModel>(interface_);

    const volScalarField CpAvg
    (
        interface_.thermo1().Cp()*interface_.thermo2().Cp()
       /(interface_.thermo1().Cp() + interface_.thermo2().Cp())
    );

    return drag.K()*CpAvg/Pr_;
}


// ************************************************************************* //
