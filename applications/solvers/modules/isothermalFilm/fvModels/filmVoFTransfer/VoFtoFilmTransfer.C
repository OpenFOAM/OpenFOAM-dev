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

#include "VoFtoFilmTransfer.H"
#include "filmToVoFTransfer.H"
#include "mappedPatchBase.H"
#include "fvmSup.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

namespace Foam
{
    namespace fv
    {
        defineTypeNameAndDebug(VoFtoFilmTransfer, 0);

        addToRunTimeSelectionTable
        (
            fvModel,
            VoFtoFilmTransfer,
            dictionary
        );
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::VoFtoFilmTransfer::VoFtoFilmTransfer
(
    const word& sourceName,
    const word& modelType,
    const fvMesh& mesh,
    const dictionary& dict
)
:
    fvModel(sourceName, modelType, mesh, dict),
    VoF_(mesh.lookupObject<solvers::compressibleVoF>(solver::typeName)),
    filmPatchName_(dict.lookup("filmPatch")),
    filmPatchi_(mesh.boundaryMesh().findPatchID(filmPatchName_)),
    phaseName_(dict.lookup("phase")),
    thermo_
    (
        phaseName_ == VoF_.mixture.phase1Name()
      ? VoF_.mixture.thermo1()
      : VoF_.mixture.thermo2()
    ),
    alpha_
    (
        phaseName_ == VoF_.mixture.phase1Name()
      ? VoF_.mixture.alpha1()
      : VoF_.mixture.alpha2()
    ),
    curTimeIndex_(-1),
    deltaFactorToFilm_
    (
        dict.lookupOrDefault<scalar>("deltaFactorToFilm", 0.5)
    ),
    alphaToFilm_
    (
        dict.lookupOrDefault<scalar>("alphaToFilm", 0.1)
    ),
    transferRateCoeff_
    (
        dict.lookupOrDefault<scalar>("transferRateCoeff", 0.1)
    ),
    transferRate_
    (
        volScalarField::Internal::New
        (
            "transferRate",
            mesh,
            dimensionedScalar(dimless/dimTime, 0)
        )
    )
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::wordList Foam::fv::VoFtoFilmTransfer::addSupFields() const
{
    return wordList
    (
        {
            alpha_.name(),
            thermo_.rho()().name(),
            thermo_.he().name(),
            VoF_.U.name()
        }
    );
}


void Foam::fv::VoFtoFilmTransfer::correct()
{
    if (curTimeIndex_ == mesh().time().timeIndex())
    {
        return;
    }

    curTimeIndex_ = mesh().time().timeIndex();


    const scalar deltaT = mesh().time().deltaTValue();

    const polyPatch& VoFFilmPatch = mesh().boundaryMesh()[filmPatchi_];


    // VoF properties

    const scalarField& alpha = alpha_.boundaryField()[filmPatchi_];

    const scalarField& deltaCoeffs =
        mesh().boundary()[filmPatchi_].deltaCoeffs();

    const labelList& faceCells = mesh().boundary()[filmPatchi_].faceCells();


    // Film properties

    const mappedPatchBase& VoFFilmPatchMap = refCast<const mappedPatchBase>
    (
        VoFFilmPatch
    );

    const solvers::isothermalFilm& film_
    (
        VoFFilmPatchMap.nbrMesh().lookupObject<solvers::isothermalFilm>
        (
            solver::typeName
        )
    );

    const label filmVoFPatchi = VoFFilmPatchMap.nbrPolyPatch().index();

    const scalarField delta
    (
        VoFFilmPatchMap.fromNeigbour
        (
            film_.delta.boundaryField()[filmVoFPatchi]
        )
    );

    transferRate_ = Zero;

    forAll(faceCells, facei)
    {
        const label celli = faceCells[facei];

        if
        (
            alpha[facei] > 0
         && delta[facei] < 2*deltaFactorToFilm_/deltaCoeffs[facei]
         && alpha[facei] < alphaToFilm_
        )
        {
            transferRate_[celli] = -transferRateCoeff_/deltaT;
        }
    }
}


template<class Type, class TransferRateFunc>
Foam::tmp<Foam::VolInternalField<Type>>
inline Foam::fv::VoFtoFilmTransfer::filmToVoFTransferRate
(
    TransferRateFunc transferRateFunc,
    const dimensionSet& dimProp
) const
{
    const mappedPatchBase& VoFFilmPatchMap = refCast<const mappedPatchBase>
    (
        mesh().boundaryMesh()[filmPatchi_]
    );

    const Foam::fvModels& fvModels
    (
        fvModels::New
        (
            refCast<const fvMesh>(VoFFilmPatchMap.nbrMesh())
        )
    );

    const filmToVoFTransfer* filmToVoFPtr = nullptr;

    forAll(fvModels, i)
    {
        if (isType<filmToVoFTransfer>(fvModels[i]))
        {
            filmToVoFPtr = &refCast<const filmToVoFTransfer>(fvModels[i]);
        }
    }

    if (!filmToVoFPtr)
    {
        FatalErrorInFunction
            << "Cannot find filmToVoFTransfer fvModel for the film region "
            << VoFFilmPatchMap.nbrMesh().name()
            << exit(FatalError);
    }

    tmp<VolInternalField<Type>> tSu
    (
        VolInternalField<Type>::New
        (
            "Su",
            mesh(),
            dimensioned<Type>(dimProp/dimTime, Zero)
        )
    );

    UIndirectList<Type>(tSu.ref(), mesh().boundary()[filmPatchi_].faceCells()) =
        VoFFilmPatchMap.fromNeigbour
        (
            (filmToVoFPtr->*transferRateFunc)()
        );

    return tSu/mesh().V();
}


void Foam::fv::VoFtoFilmTransfer::addSup
(
    fvMatrix<scalar>& eqn,
    const word& fieldName
) const
{
    if (debug)
    {
        Info<< type() << ": applying source to " << eqn.psi().name() << endl;
    }

    if (fieldName == alpha_.name())
    {
        eqn +=
            filmToVoFTransferRate<scalar>
            (
                &filmToVoFTransfer::transferRate,
                dimVolume
            )
          + fvm::Sp(transferRate_, eqn.psi());
    }
    else
    {
        FatalErrorInFunction
            << "Support for field " << fieldName << " is not implemented"
            << exit(FatalError);
    }
}


void Foam::fv::VoFtoFilmTransfer::addSup
(
    const volScalarField& alpha,
    fvMatrix<scalar>& eqn,
    const word& fieldName
) const
{
    if (debug)
    {
        Info<< type() << ": applying source to " << eqn.psi().name() << endl;
    }

    if (fieldName == thermo_.rho()().name())
    {
        eqn +=
            filmToVoFTransferRate<scalar>
            (
                &filmToVoFTransfer::rhoTransferRate,
                dimMass
            )
          + fvm::Sp(alpha()*transferRate_, eqn.psi());
    }
    else
    {
        FatalErrorInFunction
            << "Support for field " << fieldName << " is not implemented"
            << exit(FatalError);
    }
}


void Foam::fv::VoFtoFilmTransfer::addSup
(
    const volScalarField& alpha,
    const volScalarField& rho,
    fvMatrix<scalar>& eqn,
    const word& fieldName
) const
{
    if (debug)
    {
        Info<< type() << ": applying source to " << eqn.psi().name() << endl;
    }

    if (fieldName == thermo_.he().name())
    {
        eqn +=
            filmToVoFTransferRate<scalar>
            (
                &filmToVoFTransfer::heTransferRate,
                dimEnergy
            )
          + fvm::Sp(alpha()*rho()*transferRate_, eqn.psi());
    }
    else
    {
        FatalErrorInFunction
            << "Support for field " << fieldName << " is not implemented"
            << exit(FatalError);
    }
}


void Foam::fv::VoFtoFilmTransfer::addSup
(
    const volScalarField& rho,
    fvMatrix<vector>& eqn,
    const word& fieldName
) const
{
    if (debug)
    {
        Info<< type() << ": applying source to " << eqn.psi().name() << endl;
    }

    eqn +=
        filmToVoFTransferRate<vector>
        (
            &filmToVoFTransfer::UTransferRate,
            dimMass*dimVelocity
        )
      + fvm::Sp(alpha_()*thermo_.rho()()*transferRate_, eqn.psi());
}


template<class Type, class FieldType>
inline Foam::tmp<Foam::Field<Type>> Foam::fv::VoFtoFilmTransfer::TransferRate
(
    const FieldType& f
) const
{
    const labelList& faceCells = mesh().boundary()[filmPatchi_].faceCells();

    return tmp<Field<Type>>
    (
        new Field<Type>
        (
            UIndirectList<Type>
            (
                -alpha_()*transferRate_*mesh().V()*f,
                faceCells
            )
        )
    );
}


Foam::tmp<Foam::scalarField>
Foam::fv::VoFtoFilmTransfer::rhoTransferRate() const
{
    return TransferRate<scalar>(thermo_.rho()());
}


Foam::tmp<Foam::scalarField>
Foam::fv::VoFtoFilmTransfer::heTransferRate() const
{
    return TransferRate<scalar>(thermo_.rho()()*thermo_.he()());
}


Foam::tmp<Foam::vectorField>
Foam::fv::VoFtoFilmTransfer::UTransferRate() const
{
    return TransferRate<vector>(thermo_.rho()()*VoF_.U());
}


void Foam::fv::VoFtoFilmTransfer::topoChange(const polyTopoChangeMap&)
{
    transferRate_.setSize(mesh().nCells());
}


void Foam::fv::VoFtoFilmTransfer::mapMesh(const polyMeshMap& map)
{
    transferRate_.setSize(mesh().nCells());
}


void Foam::fv::VoFtoFilmTransfer::distribute(const polyDistributionMap&)
{
    transferRate_.setSize(mesh().nCells());
}


bool Foam::fv::VoFtoFilmTransfer::movePoints()
{
    return true;
}


// ************************************************************************* //
