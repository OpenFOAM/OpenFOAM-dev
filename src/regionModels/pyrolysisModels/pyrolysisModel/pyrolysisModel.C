/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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

#include "pyrolysisModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace pyrolysisModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(pyrolysisModel, 0);
defineRunTimeSelectionTable(pyrolysisModel, mesh);
defineRunTimeSelectionTable(pyrolysisModel, dictionary);

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void pyrolysisModel::readPyrolysisControls()
{}


bool pyrolysisModel::read()
{
    if (regionModel1D::read())
    {
        readPyrolysisControls();
        return true;
    }
    else
    {
        return false;
    }
}


bool pyrolysisModel::read(const dictionary& dict)
{
    if (regionModel1D::read(dict))
    {
        readPyrolysisControls();
        return true;
    }
    else
    {
        return false;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

pyrolysisModel::pyrolysisModel(const fvMesh& mesh, const word& regionType)
:
    regionModel1D(mesh, regionType)
{}


pyrolysisModel::pyrolysisModel
(
    const word& modelType,
    const fvMesh& mesh,
    const word& regionType
)
:
    regionModel1D(mesh, regionType, modelType)
{
    if (active_)
    {
        read();
    }
}


pyrolysisModel::pyrolysisModel
(
    const word& modelType,
    const fvMesh& mesh,
    const dictionary& dict,
    const word& regionType
)
:
    regionModel1D(mesh, regionType, modelType, dict)
{
    if (active_)
    {
        read(dict);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

pyrolysisModel::~pyrolysisModel()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

scalar pyrolysisModel::addMassSources
(
    const label patchi,
    const label facei
)
{
    return 0.0;
}


scalar pyrolysisModel::solidRegionDiffNo() const
{
    return -great;
}


scalar pyrolysisModel::maxDiff() const
{
    return great;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace pyrolysisModels
} // End namespace regionModels
} // End namespace Foam

// ************************************************************************* //
