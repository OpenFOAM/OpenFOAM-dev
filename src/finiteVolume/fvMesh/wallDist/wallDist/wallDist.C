/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2015-2021 OpenFOAM Foundation
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

#include "wallDist.H"
#include "wallPolyPatch.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(wallDist, 0);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::wallDist::constructn() const
{
    n_ = volVectorField::New
    (
        "n" & patchTypeName_,
        mesh(),
        dimensionedVector(dimless, Zero),
        patchDistMethod::patchTypes<vector>(mesh(), patchIDs_)
    );

    const fvPatchList& patches = mesh().boundary();

    volVectorField::Boundary& nbf = n_.ref().boundaryFieldRef();

    forAllConstIter(labelHashSet, patchIDs_, iter)
    {
        label patchi = iter.key();
        nbf[patchi] == patches[patchi].nf();
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::wallDist::wallDist(const fvMesh& mesh, const word& patchTypeName)
:
    MeshObject<fvMesh, Foam::UpdateableMeshObject, wallDist>(mesh),
    patchIDs_(mesh.boundaryMesh().findPatchIDs<wallPolyPatch>()),
    patchTypeName_(patchTypeName),
    pdm_
    (
        patchDistMethod::New
        (
            static_cast<const fvSchemes&>(mesh)
           .subDict(patchTypeName_ & "Dist"),
            mesh,
            patchIDs_
        )
    ),
    y_
    (
        IOobject
        (
            "y" & patchTypeName_,
            mesh.time().timeName(),
            mesh
        ),
        mesh,
        dimensionedScalar("y" & patchTypeName_, dimLength, small),
        patchDistMethod::patchTypes<scalar>(mesh, patchIDs_)
    ),
    nRequired_
    (
        static_cast<const fvSchemes&>(mesh).subDict(patchTypeName_ & "Dist")
       .lookupOrDefault<Switch>("nRequired", false)
    ),
    n_(volVectorField::null())
{
    if (nRequired_)
    {
        constructn();
    }

    movePoints();
}


Foam::wallDist::wallDist
(
    const fvMesh& mesh,
    const labelHashSet& patchIDs,
    const word& patchTypeName
)
:
    MeshObject<fvMesh, Foam::UpdateableMeshObject, wallDist>(mesh),
    patchIDs_(patchIDs),
    patchTypeName_(patchTypeName),
    pdm_
    (
        patchDistMethod::New
        (
            static_cast<const fvSchemes&>(mesh)
           .subDict(patchTypeName_ & "Dist"),
            mesh,
            patchIDs_
        )
    ),
    y_
    (
        IOobject
        (
            "y" & patchTypeName_,
            mesh.time().timeName(),
            mesh
        ),
        mesh,
        dimensionedScalar("y" & patchTypeName_, dimLength, small),
        patchDistMethod::patchTypes<scalar>(mesh, patchIDs_)
    ),
    nRequired_
    (
        static_cast<const fvSchemes&>(mesh).subDict(patchTypeName_ & "Dist")
       .lookupOrDefault<Switch>("nRequired", false)
    ),
    n_(volVectorField::null())
{
    if (nRequired_)
    {
        constructn();
    }

    movePoints();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::wallDist::~wallDist()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::volVectorField& Foam::wallDist::n() const
{
    if (isNull(n_()))
    {
        WarningInFunction
            << "n requested but 'nRequired' not specified in the "
            << (patchTypeName_ & "Dist") << " dictionary" << nl
            << "    Recalculating y and n fields." << endl;

        nRequired_ = true;
        constructn();
        pdm_->correct(y_, n_.ref());
    }

    return n_();
}


void Foam::wallDist::updateMesh(const mapPolyMesh& map)
{
    pdm_->updateMesh(map);
    movePoints();
}


void Foam::wallDist::distribute(const mapDistributePolyMesh& map)
{
    pdm_->distribute(map);
    movePoints();
}


bool Foam::wallDist::movePoints()
{
    if (pdm_->movePoints())
    {
        if (nRequired_)
        {
            return pdm_->correct(y_, n_.ref());
        }
        else
        {
            return pdm_->correct(y_);
        }
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
