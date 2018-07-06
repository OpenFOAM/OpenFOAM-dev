/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2018 OpenFOAM Foundation
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
#include "phasePair.H"
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
    const phasePair& pair
)
:
    phaseTransferModel(dict, pair),
    dropletName_(dict.lookup("droplet")),
    surfaceName_(dict.lookup("surface")),
    efficiency_(readScalar(dict.lookup("efficiency")))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::phaseTransferModels::deposition::~deposition()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::phaseTransferModels::deposition::dmdt() const
{
    const phaseModel* dropletPtr = nullptr;
    scalar sign = 1;
    if (dropletName_ == pair_.first())
    {
        dropletPtr = &pair_.phase1();
        sign = -1;
    }
    else if (dropletName_ == pair_.second())
    {
        dropletPtr = &pair_.phase2();
        sign = 1;
    }
    else
    {
        FatalErrorInFunction
            << "The specified droplet phase, " << dropletName_ << ", is not in "
            << "the " << pair_ << " pair"
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
