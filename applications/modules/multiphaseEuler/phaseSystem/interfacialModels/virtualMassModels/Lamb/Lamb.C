/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2014-2025 OpenFOAM Foundation
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

#include "Lamb.H"
#include "aspectRatioModel.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace virtualMassModels
{
    defineTypeNameAndDebug(Lamb, 0);
    addToRunTimeSelectionTable
    (
        virtualMassModel,
        Lamb,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::virtualMassModels::Lamb::Lamb
(
    const dictionary& dict,
    const phaseInterface& interface,
    const bool registerObject
)
:
    dispersedVirtualMassModel(dict, interface, registerObject),
    aspectRatio_(aspectRatioModel::New(dict.subDict("aspectRatio"), interface))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::virtualMassModels::Lamb::~Lamb()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::virtualMassModels::Lamb::Cvm() const
{
    volScalarField E(min(max(aspectRatio_->E(), small), 1 - small));
    volScalarField rtOmEsq(sqrt(1 - sqr(E)));

    return
        (rtOmEsq - E*acos(E))
       /(E*acos(E) - sqr(E)*rtOmEsq);
}


// ************************************************************************* //
