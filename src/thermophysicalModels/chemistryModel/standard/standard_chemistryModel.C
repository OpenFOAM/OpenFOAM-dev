/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2022-2026 OpenFOAM Foundation
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

#include "standard_chemistryModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace chemistryModels
{
    const Foam::NamedEnum<standard::jacobianType, 2>
    standard::jacobianTypeNames
    {
        "fast",
        "exact"
    };

    defineTypeNameAndDebug(standard, 0);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::chemistryModels::standard::standard
(
    const fluidMulticomponentThermo& thermo
)
:
    chemistryModel(thermo),
    ODESystem(),
    Yvf_(this->thermo().Y()),
    nSpecie_(Yvf_.size()),
    reduction_(false),
    cTos_(nSpecie_, -1),
    sToc_(nSpecie_)
{
    Info<< "chemistryModels::standard: Number of species = "
        << nSpecie_ << endl;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::chemistryModels::standard::~standard()
{}


// ************************************************************************* //
