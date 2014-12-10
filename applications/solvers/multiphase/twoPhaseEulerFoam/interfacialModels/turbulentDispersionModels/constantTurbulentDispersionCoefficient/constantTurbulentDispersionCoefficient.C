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

#include "constantTurbulentDispersionCoefficient.H"
#include "phasePair.H"
#include "fvc.H"
#include "PhaseCompressibleTurbulenceModel.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace turbulentDispersionModels
{
    defineTypeNameAndDebug(constantTurbulentDispersionCoefficient, 0);
    addToRunTimeSelectionTable
    (
        turbulentDispersionModel,
        constantTurbulentDispersionCoefficient,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::turbulentDispersionModels::constantTurbulentDispersionCoefficient::
constantTurbulentDispersionCoefficient
(
    const dictionary& dict,
    const phasePair& pair
)
:
    turbulentDispersionModel(dict, pair),
    Ctd_("Ctd", dimless, dict.lookup("Ctd"))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::turbulentDispersionModels::constantTurbulentDispersionCoefficient::
~constantTurbulentDispersionCoefficient()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volVectorField>
Foam::turbulentDispersionModels::constantTurbulentDispersionCoefficient::
F() const
{
    return
        Ctd_
       *pair_.dispersed()
       *pair_.continuous().rho()
       *pair_.continuous().turbulence().k()
       *fvc::grad(pair_.dispersed());
}


// ************************************************************************* //
