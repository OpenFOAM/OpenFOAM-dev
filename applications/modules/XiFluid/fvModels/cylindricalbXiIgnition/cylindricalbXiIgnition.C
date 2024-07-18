/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2024 OpenFOAM Foundation
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

#include "cylindricalbXiIgnition.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

namespace Foam
{
    namespace fv
    {
        defineTypeNameAndDebug(cylindricalbXiIgnition, 0);

        addToRunTimeSelectionTable
        (
            fvModel,
            cylindricalbXiIgnition,
            dictionary
        );
    }
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::fv::cylindricalbXiIgnition::readCoeffs()
{
    const dictionary& XiCorrCoeffs(coeffs().subDict("XiCorr"));

    thickness_.read(XiCorrCoeffs);
    cylinderFraction_.readIfPresent(XiCorrCoeffs);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::cylindricalbXiIgnition::cylindricalbXiIgnition
(
    const word& name,
    const word& modelType,
    const fvMesh& mesh,
    const dictionary& dict
)
:
    constantbXiIgnition(name, modelType, mesh, dict),
    thickness_("thickness", dimLength, 0),
    cylinderFraction_("cylinderFraction", dimless, 1)
{
    readCoeffs();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::dimensionedScalar Foam::fv::cylindricalbXiIgnition::Ak
(
    const dimensionedScalar& Vk
) const
{
    // Radius of the ignition kernel
    dimensionedScalar rk
    (
        sqrt(Vk/(cylinderFraction_*constant::mathematical::pi*thickness_))
    );

    // Return area of the ignition kernel
    return cylinderFraction_*2*constant::mathematical::pi*rk*thickness_;
}


bool Foam::fv::cylindricalbXiIgnition::read(const dictionary& dict)
{
    if (constantbXiIgnition::read(dict))
    {
        readCoeffs();
        return true;
    }
    else
    {
        return false;
    }

    return false;
}


// ************************************************************************* //
