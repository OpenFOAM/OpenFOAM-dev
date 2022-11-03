/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2014-2022 OpenFOAM Foundation
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

#include "Frank.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace wallLubricationModels
{
    defineTypeNameAndDebug(Frank, 0);
    addToRunTimeSelectionTable
    (
        wallLubricationModel,
        Frank,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::wallLubricationModels::Frank::Frank
(
    const dictionary& dict,
    const phaseInterface& interface
)
:
    dispersedWallLubricationModel(dict, interface),
    Cwd_("Cwd", dimless, dict),
    Cwc_("Cwc", dimless, dict),
    p_(dict.lookup<scalar>("p"))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::wallLubricationModels::Frank::~Frank()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volVectorField> Foam::wallLubricationModels::Frank::Fi() const
{
    const volVectorField Ur(interface_.Ur());

    const volVectorField& n(nWall());
    const volScalarField& y(yWall());

    const volScalarField Eo(interface_.Eo());
    const volScalarField yTilde(y/(Cwc_*interface_.dispersed().d()));

    return zeroGradWalls
    (
        (
            pos0(Eo - 1)*neg(Eo - 5)*exp(-0.933*Eo + 0.179)
          + pos0(Eo - 5)*neg(Eo - 33)*(0.00599*Eo - 0.0187)
          + pos0(Eo - 33)*0.179
        )
       *max
        (
            dimensionedScalar(dimless/dimLength, 0),
            (1 - yTilde)/(Cwd_*y*pow(yTilde, p_ - 1))
        )
       *interface_.continuous().rho()
       *magSqr(Ur - (Ur & n)*n)
       *n
    );
}


// ************************************************************************* //
