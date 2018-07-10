/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2017-2018 OpenFOAM Foundation
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

#include "rhoThermo.H"
#include "heRhoThermo.H"
#include "pureMixture.H"
#include "thermo.H"
#include "sensibleInternalEnergy.H"
#include "sensibleEnthalpy.H"
#include "thermophysicalPropertiesSelector.H"
#include "liquidProperties.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/* * * * * * * * * * * * * * * private static data * * * * * * * * * * * * * */

typedef heRhoThermo
<
    rhoThermo,
    pureMixture
    <
        species::thermo
        <
            thermophysicalPropertiesSelector<liquidProperties>,
            sensibleInternalEnergy
        >
    >
> heRhoThermopureMixtureliquidProperties;


typedef heRhoThermo
<
    rhoThermo,
    pureMixture
    <
        species::thermo
        <
            thermophysicalPropertiesSelector<liquidProperties>,
            sensibleEnthalpy
        >
    >
> heRhoThermopureMixtureEnthalpyliquidProperties;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
