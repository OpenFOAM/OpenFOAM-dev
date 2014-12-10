/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2014 OpenFOAM Foundation
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

#include "constantLiftCoefficient.H"
#include "phasePair.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace liftModels
{
    defineTypeNameAndDebug(constantLiftCoefficient, 0);
    addToRunTimeSelectionTable(liftModel, constantLiftCoefficient, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::liftModels::constantLiftCoefficient::constantLiftCoefficient
(
    const dictionary& dict,
    const phasePair& pair
)
:
    liftModel(dict, pair),
    Cl_("Cl", dimless, dict.lookup("Cl"))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::liftModels::constantLiftCoefficient::~constantLiftCoefficient()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::liftModels::constantLiftCoefficient::Cl() const
{
    const fvMesh& mesh(this->pair_.phase1().mesh());

    return
        tmp<volScalarField>
        (
            new volScalarField
            (
                IOobject
                (
                    "zero",
                    mesh.time().timeName(),
                    mesh
                ),
                mesh,
                Cl_
            )
        );
}


// ************************************************************************* //
