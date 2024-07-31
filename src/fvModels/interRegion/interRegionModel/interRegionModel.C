/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2024 OpenFOAM Foundation
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
#include "fvModels.H"
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
    master_ = dict.lookupOrDefault<bool>("master", true);

    nbrRegionName_ =
        dict.lookupBackwardsCompatible<word>
        ({
            "nbrRegion",
            "nbrRegionName"
        });

    dict.lookup("interpolationMethod") >> interpolationMethod_;
}


const Foam::cellsToCells& Foam::fv::interRegionModel::interpolation() const
{
    if (!interpolationPtr_.valid())
    {
        Info<< incrIndent;

        if (master_)
        {
            Info<< indent << "- selecting inter region mapping" << endl;

            if (mesh().name() == nbrMesh().name())
            {
                FatalErrorInFunction
                    << "Inter-region model selected, but local and "
                    << "neighbour regions are the same: " << nl
                    << "    local region: " << mesh().name() << nl
                    << "    secondary region: " << nbrMesh().name() << nl
                    << exit(FatalError);
            }

            if (mesh().bounds().overlaps(nbrMesh().bounds()))
            {
                interpolationPtr_ = cellsToCells::New(interpolationMethod_);
                interpolationPtr_->update(mesh(), nbrMesh());
            }
            else
            {
                FatalErrorInFunction
                    << "regions " << mesh().name() << " and "
                    << nbrMesh().name() <<  " do not intersect"
                    << exit(FatalError);
            }
        }

        Info<< decrIndent;
    }

    return interpolationPtr_();
}


const Foam::fv::interRegionModel& Foam::fv::interRegionModel::nbrModel() const
{
    const fvMesh& nbrMesh = mesh().time().lookupObject<fvMesh>(nbrRegionName());

    const PtrListDictionary<fvModel>& fvModels =
        nbrMesh.lookupObject<Foam::fvModels>("fvModels");

    forAll(fvModels, fvModeli)
    {
        if (isA<interRegionModel>(fvModels[fvModeli]))
        {
            const interRegionModel& model =
                refCast<const interRegionModel>(fvModels[fvModeli]);

            if (model.nbrRegionName() == mesh().name())
            {
                return model;
            }
        }
    }

    FatalErrorInFunction
        << "Neighbour model not found in region " << nbrMesh.name() << nl
        << exit(FatalError);
    return NullObjectRef<interRegionModel>();
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
    master_(false),
    nbrRegionName_(word::null),
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
