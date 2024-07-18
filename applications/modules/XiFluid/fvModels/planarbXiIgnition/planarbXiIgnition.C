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

#include "planarbXiIgnition.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

namespace Foam
{
    namespace fv
    {
        defineTypeNameAndDebug(planarbXiIgnition, 0);

        addToRunTimeSelectionTable
        (
            fvModel,
            planarbXiIgnition,
            dictionary
        );
    }
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::fv::planarbXiIgnition::readCoeffs()
{
    const dictionary& XiCorrCoeffs(coeffs().subDict("XiCorr"));

    area_.read(XiCorrCoeffs);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::planarbXiIgnition::planarbXiIgnition
(
    const word& name,
    const word& modelType,
    const fvMesh& mesh,
    const dictionary& dict
)
:
    constantbXiIgnition(name, modelType, mesh, dict),
    area_("area", dimLength, 0)
{
    readCoeffs();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::dimensionedScalar Foam::fv::planarbXiIgnition::Ak
(
    const dimensionedScalar& Vk
) const
{
    return area_;
}


bool Foam::fv::planarbXiIgnition::read(const dictionary& dict)
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
