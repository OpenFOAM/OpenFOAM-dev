/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2023 OpenFOAM Foundation
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

#include "surfaceFilm.H"
#include "timeIOdictionary.H"
#include "mappedWallPolyPatch.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(surfaceFilm, 0);
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

bool Foam::surfaceFilm::read()
{
    if (regIOobject::read())
    {
        if (const dictionary* dictPtr = subDictPtr(modelType_ + "Coeffs"))
        {
            coeffs_ <<= *dictPtr;
        }

        infoOutput_.readIfPresent("infoOutput", *this);

        return true;
    }
    else
    {
        return false;
    }
}


Foam::label Foam::surfaceFilm::nbrCoupledPatchID
(
    const surfaceFilm& nbrFilm,
    const label filmPatchi
) const
{
    label nbrPatchi = -1;

    // film
    const fvMesh& nbrFilmMesh = nbrFilm.mesh();

    // boundary meshes
    const polyBoundaryMesh& pbm = mesh().boundaryMesh();
    const polyBoundaryMesh& nbrPbm = nbrFilmMesh.boundaryMesh();

    if (filmPatchi > pbm.size() - 1)
    {
        FatalErrorInFunction
            << "film patch index out of bounds: "
            << "film patch index = " << filmPatchi
            << ", maximum index = " << pbm.size() - 1
            << abort(FatalError);
    }

    const mappedPatchBase& mpb =
        mappedPatchBase::getMap(pbm[filmPatchi]);

    // sample patch name on the primary film
    const word& primaryPatchName = mpb.nbrPatchName();

    // find patch on nbr film that has the same sample patch name
    forAll(nbrFilm.intCoupledPatchIDs(), j)
    {
        const label nbrFilmPatchi = nbrFilm.intCoupledPatchIDs()[j];

        const mappedPatchBase& mpb =
            mappedPatchBase::getMap(nbrPbm[nbrFilmPatchi]);

        if (mpb.nbrPatchName() == primaryPatchName)
        {
            nbrPatchi = nbrFilmPatchi;
            break;
        }
    }

    if (nbrPatchi == -1)
    {
        const polyPatch& p = mesh().boundaryMesh()[filmPatchi];

        FatalErrorInFunction
            << "Unable to find patch pair for local patch "
            << p.name() << " and film " << nbrFilm.name()
            << abort(FatalError);
    }

    return nbrPatchi;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::surfaceFilm::surfaceFilm
(
    const word& modelType,
    const fvMesh& primaryMesh,
    const dimensionedVector& g,
    const word& filmType
)
:
    IOdictionary
    (
        IOobject
        (
            filmType + "Properties",
            primaryMesh.time().constant(),
            primaryMesh.time(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    primaryMesh_(primaryMesh),
    time_(primaryMesh.time()),
    infoOutput_(true),
    modelType_(modelType),
    regionName_(lookup("regionName")),
    mesh_
    (
        IOobject
        (
            regionName_,
            time_.name(),
            time_,
            IOobject::MUST_READ
        )
    ),
    coeffs_(optionalSubDict(modelType + "Coeffs")),
    outputPropertiesPtr_(nullptr),
    primaryPatchIDs_(),
    intCoupledPatchIDs_(),
    passivePatchIDs_(),
    nHat_
    (
        IOobject
        (
            "nHat",
            time_.name(),
            mesh_
        ),
        mesh_,
        dimensionedVector(dimless, Zero),
        zeroGradientFvPatchField<vector>::typeName
    ),
    magSf_
    (
        IOobject
        (
            "magSf",
            time_.name(),
            mesh_
        ),
        mesh_,
        dimensionedScalar(dimArea, 0)
    ),
    VbyA_
    (
        IOobject
        (
            "VbyA",
            time_.name(),
            mesh_
        ),
        mesh_,
        dimensionedScalar(dimLength, 0),
        zeroGradientFvPatchField<vector>::typeName
    ),
    g_(g)
{
    label nBoundaryFaces = 0;
    DynamicList<label> primaryPatchIDs;
    DynamicList<label> intCoupledPatchIDs;
    const polyBoundaryMesh& rbm = mesh_.boundaryMesh();

    forAll(rbm, patchi)
    {
        const polyPatch& filmPatch = rbm[patchi];

        if (isA<mappedPatchBase>(filmPatch))
        {
            intCoupledPatchIDs.append(patchi);

            nBoundaryFaces += filmPatch.faceCells().size();

            const mappedPatchBase& mapPatch =
                refCast<const mappedPatchBase>(filmPatch);

            if
            (
                primaryMesh_.time().foundObject<polyMesh>
                (
                    mapPatch.nbrRegionName()
                )
            )
            {
                const label primaryPatchi = mapPatch.nbrPolyPatch().index();
                primaryPatchIDs.append(primaryPatchi);
            }
        }
    }

    primaryPatchIDs_.transfer(primaryPatchIDs);
    intCoupledPatchIDs_.transfer(intCoupledPatchIDs);

    if (returnReduce(nBoundaryFaces, sumOp<label>()) == 0)
    {
        WarningInFunction
            << "Film model has no mapped boundary conditions - transfer "
            << "between films will not be possible" << endl;
    }

    if (!outputPropertiesPtr_.valid())
    {
        outputPropertiesPtr_.reset
        (
            new timeIOdictionary
            (
                IOobject
                (
                    regionName_ + "OutputProperties",
                    time_.name(),
                    "uniform"/regionName_,
                    primaryMesh_,
                    IOobject::READ_IF_PRESENT,
                    IOobject::NO_WRITE
                )
            )
        );
    }

    read();

    nBoundaryFaces = 0;

    forAll(intCoupledPatchIDs_, i)
    {
        const label patchi = intCoupledPatchIDs_[i];
        const polyPatch& pp = rbm[patchi];
        const labelList& fCells = pp.faceCells();

        nBoundaryFaces += fCells.size();

        UIndirectList<vector>(nHat_, fCells) = pp.faceNormals();
        UIndirectList<scalar>(magSf_, fCells) = pp.magFaceAreas();
    }
    nHat_.correctBoundaryConditions();

    if (nBoundaryFaces != mesh_.nCells())
    {
        FatalErrorInFunction
            << "Number of primary film coupled boundary faces not equal to "
            << "the number of cells in the local film" << nl << nl
            << "Number of cells = " << mesh_.nCells() << nl
            << "Boundary faces  = " << nBoundaryFaces << nl
            << abort(FatalError);
    }

    passivePatchIDs_.setSize(intCoupledPatchIDs_.size(), -1);
    forAll(intCoupledPatchIDs_, i)
    {
        const label patchi = intCoupledPatchIDs_[i];
        const polyPatch& ppIntCoupled = rbm[patchi];
        if (ppIntCoupled.size() > 0)
        {
            const label cellId = rbm[patchi].faceCells()[0];
            const cell& cFaces = mesh_.cells()[cellId];
            const label faceO
            (
                cFaces.opposingFaceLabel
                (
                    ppIntCoupled.start(), mesh_.faces()
                )
            );
            passivePatchIDs_[i] = rbm.whichPatch(faceO);
        }
    }

    Pstream::listCombineGather(passivePatchIDs_, maxEqOp<label>());
    Pstream::listCombineScatter(passivePatchIDs_);

    VbyA_.primitiveFieldRef() = mesh_.V()/magSf_;
    VbyA_.correctBoundaryConditions();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::surfaceFilm::~surfaceFilm()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

const Foam::volVectorField&
Foam::surfaceFilm::nHat() const
{
    return nHat_;
}


const Foam::volScalarField::Internal&
Foam::surfaceFilm::magSf() const
{
    return magSf_;
}


const Foam::volScalarField&
Foam::surfaceFilm::VbyA() const
{
    return VbyA_;
}


const Foam::labelList&
Foam::surfaceFilm::passivePatchIDs() const
{
    return passivePatchIDs_;
}


void Foam::surfaceFilm::evolve()
{
    Info<< nl << "Evolving " << modelType_ << " for film "
        << mesh_.name() << endl;

    preEvolveFilm();

    evolveFilm();

    postEvolveFilm();

    // Provide some feedback
    if (infoOutput_)
    {
        Info<< incrIndent;
        info();
        Info<< endl << decrIndent;
    }

    if (time_.writeTime())
    {
        outputProperties().writeObject
        (
            IOstream::ASCII,
            IOstream::currentVersion,
            time_.writeCompression(),
            true
        );
    }
}


void Foam::surfaceFilm::preEvolveFilm()
{}


void Foam::surfaceFilm::evolveFilm()
{}


void Foam::surfaceFilm::postEvolveFilm()
{}


// ************************************************************************* //
