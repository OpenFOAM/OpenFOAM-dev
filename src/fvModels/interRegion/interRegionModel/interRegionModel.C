/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2026 OpenFOAM Foundation
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

#include "interRegionModel.H"
#include "matchingCellsToCells.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(interRegionModel, 0);
}
}


// * * * * * * * * * * * * *  Private Member Functions * * * * * * * * * * * //

void Foam::fv::interRegionModel::readCoeffs(const dictionary& dict)
{
    owner_ =
        dict.found("owner") || dict.found("master")
      ? dict.lookupBackwardsCompatible<bool>({"owner", "master"})
      : -1;

    nbrRegionName_ =
        dict.lookupBackwardsCompatible<word>
        ({
            "neighbourRegion",
            "nbrRegion",
            "nbrRegionName"
        });

    if (nbrRegionName_ == mesh().name())
    {
        FatalIOErrorInFunction(dict)
            << "Neighbour region is the same as the region"
            << exit(FatalError);
    }

    nbrModelName_ =
        dict.lookupOrDefaultBackwardsCompatible<word>
        (
            {"neighbourModel", "nbrModel"},
            word::null
        );

    nbrModelNameDict_ = dict;

    dict.lookup("interpolationMethod") >> interpolationMethod_;
}


const Foam::cellsToCells& Foam::fv::interRegionModel::interpolation() const
{
    if (neighbour())
    {
        FatalErrorInFunction
            << "Inter-region mapping is not available to the neighbour model"
            << exit(FatalError);
    }

    if (!interpolationPtr_.valid())
    {
        Info<< incrIndent;

        Info<< indent << "- selecting inter region mapping" << endl;

        interpolationPtr_ = cellsToCells::New(interpolationMethod_);
        interpolationPtr_->update(mesh(), nbrMesh());

        Info<< decrIndent;
    }

    return interpolationPtr_();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::interRegionModel::interRegionModel
(
    const word& name,
    const word& modelType,
    const fvMesh& mesh,
    const dictionary& dict
)
:
    fvModel(name, modelType, mesh, dict),
    owner_(-1),
    nbrRegionName_(word::null),
    nbrModelName_(word::null),
    nbrModelNameDict_(),
    interpolationMethod_(cellsToCellss::matching::typeName),
    interpolationPtr_()
{
    readCoeffs(coeffs(dict));
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fv::interRegionModel::~interRegionModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::fv::interRegionModel::read(const dictionary& dict)
{
    if (fvModel::read(dict))
    {
        readCoeffs(coeffs(dict));
        interpolationPtr_.clear();
        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
