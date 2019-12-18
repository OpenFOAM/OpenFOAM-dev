/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
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

#include "thermocapillaryForce.H"
#include "addToRunTimeSelectionTable.H"
#include "fvcGrad.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace surfaceFilmModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(thermocapillaryForce, 0);
addToRunTimeSelectionTable(force, thermocapillaryForce, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

thermocapillaryForce::thermocapillaryForce
(
    surfaceFilmRegionModel& film,
    const dictionary& dict
)
:
    force(film)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

thermocapillaryForce::~thermocapillaryForce()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

tmp<fvVectorMatrix> thermocapillaryForce::correct(volVectorField& U)
{
    const volScalarField& sigma = filmModel_.sigma();

    tmp<fvVectorMatrix>
        tfvm(new fvVectorMatrix(U, dimForce));

    tfvm.ref() += fvc::grad(sigma)/filmModel_.VbyA();

    return tfvm;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace surfaceFilmModels
} // End namespace regionModels
} // End namespace Foam

// ************************************************************************* //
