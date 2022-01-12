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

#include "Antal.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace wallLubricationModels
{
    defineTypeNameAndDebug(Antal, 0);
    addToRunTimeSelectionTable(wallLubricationModel, Antal, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::wallLubricationModels::Antal::Antal
(
    const dictionary& dict,
    const phaseInterface& interface
)
:
    dispersedWallLubricationModel(dict, interface),
    Cw1_("Cw1", dimless, dict),
    Cw2_("Cw2", dimless, dict)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::wallLubricationModels::Antal::~Antal()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volVectorField> Foam::wallLubricationModels::Antal::Fi() const
{
    const volVectorField Ur(interface_.Ur());

    const volVectorField& n(nWall());

    return zeroGradWalls
    (
        max
        (
            dimensionedScalar(dimless/dimLength, 0),
            Cw1_/interface_.dispersed().d() + Cw2_/yWall()
        )
       *interface_.continuous().rho()
       *magSqr(Ur - (Ur & n)*n)
       *n
    );
}


// ************************************************************************* //
