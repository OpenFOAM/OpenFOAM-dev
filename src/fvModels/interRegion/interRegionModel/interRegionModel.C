/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2021 OpenFOAM Foundation
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

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(interRegionModel, 0);
}
}


// * * * * * * * * * * * * *  Private Member Functions * * * * * * * * * * * //

void Foam::fv::interRegionModel::readCoeffs()
{
    master_ = coeffs().lookupOrDefault<bool>("master", true);

    nbrRegionName_ =
        coeffs().lookupBackwardsCompatible<word>
        ({
            "nbrRegion",
            "nbrRegionName"
        });

    interpolationMethod_ =
        meshToMesh::interpolationMethodNames_.read
        (
            coeffs().lookup("interpolationMethod")
        );
}


void Foam::fv::interRegionModel::setMapper() const
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
            meshInterpPtr_.reset
            (
                new meshToMesh
                (
                    mesh(),
                    nbrMesh(),
                    interpolationMethod_,
                    false // not interpolating patches
                )
            );
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
    const dictionary& dict,
    const fvMesh& mesh
)
:
    fvModel(name, modelType, dict, mesh),
    master_(false),
    nbrRegionName_(word::null),
    interpolationMethod_(meshToMesh::imDirect),
    meshInterpPtr_()
{
    readCoeffs();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fv::interRegionModel::~interRegionModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::fv::interRegionModel::read(const dictionary& dict)
{
    if (fvModel::read(dict))
    {
        readCoeffs();
        setMapper();
        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
