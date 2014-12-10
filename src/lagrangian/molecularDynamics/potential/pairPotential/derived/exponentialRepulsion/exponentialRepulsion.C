/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

#include "exponentialRepulsion.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace pairPotentials
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(exponentialRepulsion, 0);

addToRunTimeSelectionTable
(
    pairPotential,
    exponentialRepulsion,
    dictionary
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

exponentialRepulsion::exponentialRepulsion
(
    const word& name,
    const dictionary& exponentialRepulsion
)
:
    pairPotential(name, exponentialRepulsion),
    exponentialRepulsionCoeffs_
    (
        exponentialRepulsion.subDict(typeName + "Coeffs")
    ),
    rm_(readScalar(exponentialRepulsionCoeffs_.lookup("rm"))),
    epsilon_(readScalar(exponentialRepulsionCoeffs_.lookup("epsilon")))
{
    setLookupTables();
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

scalar exponentialRepulsion::unscaledEnergy(const scalar r) const
{
    return epsilon_ * exp(-r/rm_);
}


bool exponentialRepulsion::read(const dictionary& exponentialRepulsion)
{
    pairPotential::read(exponentialRepulsion);

    exponentialRepulsionCoeffs_ =
        exponentialRepulsion.subDict(typeName + "Coeffs");

    exponentialRepulsionCoeffs_.lookup("rm") >> rm_;
    exponentialRepulsionCoeffs_.lookup("epsilon") >> epsilon_;

    return true;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace pairPotentials
} // End namespace Foam

// ************************************************************************* //
