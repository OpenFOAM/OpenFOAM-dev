/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2023 OpenFOAM Foundation
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

#include "filmCloudTransfer.H"
#include "mappedPatchBase.H"
#include "fvmSup.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

namespace Foam
{
    namespace fv
    {
        defineTypeNameAndDebug(filmCloudTransfer, 0);

        addToRunTimeSelectionTable
        (
            fvModel,
            filmCloudTransfer,
            dictionary
        );
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::filmCloudTransfer::filmCloudTransfer
(
    const word& sourceName,
    const word& modelType,
    const fvMesh& mesh,
    const dictionary& dict
)
:
    fvModel(sourceName, modelType, mesh, dict),
    film_(mesh.lookupObject<solvers::isothermalFilm>(solver::typeName)),
    cloudFieldsTransferred_(false),
    correctEjection_(false),
    ejection_
    (
        dict.found("ejection")
      ? ejectionModel::New(dict.subDict("ejection"), film_)
      : autoPtr<ejectionModel>(nullptr)
    )
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::wordList Foam::fv::filmCloudTransfer::addSupFields() const
{
    return wordList
    {
        film_.alpha.name(),
        film_.thermo.he().name(),
        film_.U.name()
    };
}


void Foam::fv::filmCloudTransfer::correct()
{
    if (ejection_.valid() && correctEjection_)
    {
        ejection_->correct();

        // Do not correct ejection rate until the cloud has evolved
        // to include the last set of ejected parcels
        correctEjection_ = false;
    }
}


template<class Type>
Foam::tmp<Foam::VolInternalField<Type>>
inline Foam::fv::filmCloudTransfer::CloudToFilmTransferRate
(
    const Field<Type>& prop,
    const dimensionSet& dimProp
) const
{
    tmp<VolInternalField<Type>> tSu
    (
        VolInternalField<Type>::New
        (
            "Su",
            mesh(),
            dimensioned<Type>(dimProp/dimVolume/dimTime, Zero)
        )
    );

    if (cloudFieldsTransferred_)
    {
        const fvMesh& cloudMesh =
            refCast<const fvMesh>(film_.surfacePatchMap().nbrMesh());

        const label cloudPatchi =
            film_.surfacePatchMap().nbrPolyPatch().index();

        UIndirectList<Type>(tSu.ref(), film_.surfacePatch().faceCells()) =
            film_.surfacePatchMap().fromNeighbour
            (
                prop/cloudMesh.boundary()[cloudPatchi].magSf()
            );

        tSu.ref().field() /= film_.VbyA;
        tSu.ref().field() /= mesh().time().deltaTValue();
    }

    return tSu;
}


void Foam::fv::filmCloudTransfer::addSup
(
    const volScalarField& rho,
    const volScalarField& alpha,
    fvMatrix<scalar>& eqn
) const
{
    if (debug)
    {
        Info<< type() << ": applying source to " << eqn.psi().name() << endl;
    }

    if (&alpha == &film_.alpha && &eqn.psi() == &film_.alpha)
    {
        eqn += CloudToFilmTransferRate<scalar>(massFromCloud_, dimMass);

        if (ejection_.valid())
        {
            eqn -= fvm::Sp(ejection_->rate()*rho(), eqn.psi());
        }
    }
    else
    {
        FatalErrorInFunction
            << "Support for field " << alpha.name() << " is not implemented"
            << exit(FatalError);
    }
}


void Foam::fv::filmCloudTransfer::addSup
(
    const volScalarField& alpha,
    const volScalarField& rho,
    const volScalarField& he,
    fvMatrix<scalar>& eqn
) const
{
    if (debug)
    {
        Info<< type() << ": applying source to " << eqn.psi().name() << endl;
    }

    if (&he == &film_.thermo.he() && &eqn.psi() == &film_.thermo.he())
    {
        eqn += CloudToFilmTransferRate<scalar>(energyFromCloud_, dimEnergy);

        if (ejection_.valid())
        {
            eqn -= fvm::Sp(alpha()*rho()*ejection_->rate(), eqn.psi());
        }
    }
    else
    {
        FatalErrorInFunction
            << "Support for field " << he.name() << " is not implemented"
            << exit(FatalError);
    }
}


void Foam::fv::filmCloudTransfer::addSup
(
    const volScalarField& alpha,
    const volScalarField& rho,
    const volVectorField& U,
    fvMatrix<vector>& eqn
) const
{
    if (debug)
    {
        Info<< type() << ": applying source to " << eqn.psi().name() << endl;
    }

    if (&U == &film_.U && &U == &film_.U)
    {
        eqn += CloudToFilmTransferRate<vector>(momentumFromCloud_, dimMomentum);

        if (ejection_.valid())
        {
            eqn -= fvm::Sp(alpha()*rho()*ejection_->rate(), eqn.psi());
        }
    }
    else
    {
        FatalErrorInFunction
            << "Support for field " << U.name() << " is not implemented"
            << exit(FatalError);
    }
}


void Foam::fv::filmCloudTransfer::resetFromCloudFields()
{
    const fvMesh& cloudMesh =
        refCast<const fvMesh>(film_.surfacePatchMap().nbrMesh());
    const label cloudPatchi = film_.surfacePatchMap().nbrPolyPatch().index();
    const label nCloudPatchFaces = cloudMesh.boundary()[cloudPatchi].size();

    if (massFromCloud_.size() != nCloudPatchFaces)
    {
        massFromCloud_.setSize(nCloudPatchFaces);
        momentumFromCloud_.setSize(nCloudPatchFaces);
        energyFromCloud_.setSize(nCloudPatchFaces);
    }

    massFromCloud_ = 0;
    momentumFromCloud_ = Zero;
    energyFromCloud_ = 0;

    cloudFieldsTransferred_ = true;

    // Enable ejection correction on next call to correct()
    correctEjection_ = true;
}


void Foam::fv::filmCloudTransfer::parcelFromCloud
(
    const label facei,
    const scalar mass,
    const vector& momentum,
    const scalar energy
)
{
    massFromCloud_[facei] += mass;
    momentumFromCloud_[facei] += momentum;
    energyFromCloud_[facei] += energy;
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
inline Foam::fv::filmCloudTransfer::filmToCloudTransfer
(
    const VolInternalField<Type>& prop
) const
{
    return film_.surfacePatchMap().toNeighbour
    (
        Field<Type>(prop, film_.surfacePatch().faceCells())
    );
}


bool Foam::fv::filmCloudTransfer::ejecting() const
{
    return ejection_.valid();
}


Foam::tmp<Foam::Field<Foam::scalar>>
Foam::fv::filmCloudTransfer::ejectedMassToCloud() const
{
    return filmToCloudTransfer<scalar>
    (
        (
            mesh().V()
           *mesh().time().deltaTValue()
           *film_.alpha()*film_.rho()*ejection_->rate()
        )()
    );
}


Foam::tmp<Foam::Field<Foam::scalar>>
Foam::fv::filmCloudTransfer::ejectedDiameterToCloud() const
{
    return filmToCloudTransfer<scalar>(ejection_->diameter());
}


Foam::tmp<Foam::Field<Foam::scalar>>
Foam::fv::filmCloudTransfer::deltaToCloud() const
{
    return filmToCloudTransfer<scalar>(film_.delta);
}


Foam::tmp<Foam::Field<Foam::vector>>
Foam::fv::filmCloudTransfer::UToCloud() const
{
    return filmToCloudTransfer<vector>(film_.U);
}


Foam::tmp<Foam::Field<Foam::scalar>>
Foam::fv::filmCloudTransfer::rhoToCloud() const
{
    return filmToCloudTransfer<scalar>(film_.rho);
}


Foam::tmp<Foam::Field<Foam::scalar>>
Foam::fv::filmCloudTransfer::TToCloud() const
{
    return filmToCloudTransfer<scalar>(film_.thermo.T());
}


Foam::tmp<Foam::Field<Foam::scalar>>
Foam::fv::filmCloudTransfer::CpToCloud() const
{
    return filmToCloudTransfer<scalar>(film_.thermo.Cp());
}


void Foam::fv::filmCloudTransfer::topoChange(const polyTopoChangeMap& map)
{
    // Set the cloud field state to false, will be updated by the cloud tracking
    // If the film is evaluated before the cloud it would be better
    // if the cloud fields were mapped
    cloudFieldsTransferred_ = false;

    if (ejection_.valid())
    {
        ejection_->topoChange(map);
    }
}


void Foam::fv::filmCloudTransfer::mapMesh(const polyMeshMap& map)
{
    // Set the cloud field state to false, will be updated by the cloud tracking
    // If the film is evaluated before the cloud it would be better
    // if the cloud fields were mapped
    cloudFieldsTransferred_ = false;

    if (ejection_.valid())
    {
        ejection_->mapMesh(map);
    }
}


void Foam::fv::filmCloudTransfer::distribute(const polyDistributionMap& map)
{
    // Set the cloud field state to false, will be updated by the cloud tracking
    // If the film is evaluated before the cloud it would be better
    // if the cloud fields were mapped
    cloudFieldsTransferred_ = false;

    if (ejection_.valid())
    {
        ejection_->distribute(map);
    }
}


bool Foam::fv::filmCloudTransfer::movePoints()
{
    if (ejection_.valid())
    {
        ejection_->movePoints();
    }

    return true;
}


// ************************************************************************* //
