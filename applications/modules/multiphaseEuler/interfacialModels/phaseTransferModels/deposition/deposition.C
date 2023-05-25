/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2018-2023 OpenFOAM Foundation
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

#include "deposition.H"
#include "phaseSystem.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace phaseTransferModels
{
    defineTypeNameAndDebug(deposition, 0);
    addToRunTimeSelectionTable(phaseTransferModel, deposition, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::phaseTransferModels::deposition::deposition
(
    const dictionary& dict,
    const phaseInterface& interface
)
:
    phaseTransferModel(dict, interface),
    interface_(interface),
    dropletName_(dict.lookup("droplet")),
    surfaceName_(dict.lookup("surface")),
    efficiency_(dict.lookup<scalar>("efficiency"))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::phaseTransferModels::deposition::~deposition()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::phaseTransferModels::deposition::mixture() const
{
    return true;
}


Foam::tmp<Foam::volScalarField>
Foam::phaseTransferModels::deposition::dmdtf() const
{
    const phaseModel* dropletPtr = nullptr;
    scalar sign = 1;
    if (dropletName_ == interface_.phase1().name())
    {
        dropletPtr = &interface_.phase1();
        sign = -1;
    }
    else if (dropletName_ == interface_.phase2().name())
    {
        dropletPtr = &interface_.phase2();
        sign = 1;
    }
    else
    {
        FatalErrorInFunction
            << "The specified droplet phase, " << dropletName_ << ", is not in "
            << "the " << interface_ << " pair"
            << exit(FatalError);
    }

    const phaseModel& droplet = *dropletPtr;
    const phaseModel& surface = droplet.fluid().phases()[surfaceName_];

    return
        1.5*sign*efficiency_
       *droplet.rho()*droplet*surface/surface.d()
       *mag(droplet.U() - surface.U());
}


// ************************************************************************* //
