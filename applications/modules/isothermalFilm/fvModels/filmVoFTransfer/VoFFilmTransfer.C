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

#include "VoFFilmTransfer.H"
#include "filmVoFTransfer.H"
#include "mappedFvPatchBaseBase.H"
#include "fvmSup.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

namespace Foam
{
    namespace fv
    {
        defineTypeNameAndDebug(VoFFilmTransfer, 0);

        addToRunTimeSelectionTable
        (
            fvModel,
            VoFFilmTransfer,
            dictionary
        );
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::VoFFilmTransfer::VoFFilmTransfer
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
    filmPatchi_(mesh.boundaryMesh().findIndex(filmPatchName_)),
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

Foam::wordList Foam::fv::VoFFilmTransfer::addSupFields() const
{
    return wordList
    {
        thermo_.rho()().name(),
        thermo_.he().name(),
        VoF_.U.name()
    };
}


void Foam::fv::VoFFilmTransfer::correct()
{
    if (curTimeIndex_ == mesh().time().timeIndex())
    {
        return;
    }

    curTimeIndex_ = mesh().time().timeIndex();


    const scalar deltaT = mesh().time().deltaTValue();

    const fvPatch& VoFFilmPatch = mesh().boundary()[filmPatchi_];

    const scalarField& deltaCoeffs = VoFFilmPatch.deltaCoeffs();

    const labelList& faceCells = VoFFilmPatch.faceCells();


    // VoF properties

    const scalarField& alpha = alpha_.boundaryField()[filmPatchi_];


    // Film properties

    const mappedFvPatchBaseBase& VoFFilmPatchMap =
        refCast<const mappedFvPatchBaseBase>(VoFFilmPatch);

    const solvers::isothermalFilm& film_
    (
        VoFFilmPatchMap.nbrMesh().lookupObject<solvers::isothermalFilm>
        (
            solver::typeName
        )
    );

    const label filmVoFPatchi = VoFFilmPatchMap.nbrFvPatch().index();

    const scalarField delta
    (
        VoFFilmPatchMap.fromNeighbour
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
            transferRate_[celli] = transferRateCoeff_/deltaT;
        }
    }
}


template<class Type, class TransferRateFunc>
Foam::tmp<Foam::VolInternalField<Type>>
inline Foam::fv::VoFFilmTransfer::filmVoFTransferRate
(
    TransferRateFunc transferRateFunc,
    const dimensionSet& dimProp
) const
{
    const fvPatch& VoFFilmPatch = mesh().boundary()[filmPatchi_];

    const labelList& faceCells = VoFFilmPatch.faceCells();


    const mappedFvPatchBaseBase& VoFFilmPatchMap =
        refCast<const mappedFvPatchBaseBase>(VoFFilmPatch);

    const Foam::fvModels& fvModels =
        fvModels::New(VoFFilmPatchMap.nbrMesh());

    const filmVoFTransfer* filmVoFPtr = nullptr;

    forAll(fvModels, i)
    {
        if (isType<filmVoFTransfer>(fvModels[i]))
        {
            filmVoFPtr = &refCast<const filmVoFTransfer>(fvModels[i]);
        }
    }

    if (!filmVoFPtr)
    {
        FatalErrorInFunction
            << "Cannot find filmVoFTransfer fvModel for the film region "
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

    UIndirectList<Type>(tSu.ref(), faceCells) =
        VoFFilmPatchMap.fromNeighbour
        (
            (filmVoFPtr->*transferRateFunc)()
        );

    return tSu/mesh().V();
}


void Foam::fv::VoFFilmTransfer::addSup
(
    const volScalarField& alpha,
    const volScalarField& rho,
    fvMatrix<scalar>& eqn
) const
{
    if (debug)
    {
        Info<< type() << ": applying source to " << eqn.psi().name() << endl;
    }

    if (&rho == &thermo_.rho()())
    {
        // Explicit transfer into the VoF
        eqn +=
            filmVoFTransferRate<scalar>
            (
                &filmVoFTransfer::rhoTransferRate,
                dimMass
            );

        // Potentially implicit transfer out of the VoF
        if (&rho == &eqn.psi())
        {
            eqn -= fvm::Sp(alpha()*transferRate_, eqn.psi());
        }
        else if (&alpha == &eqn.psi())
        {
            eqn -= fvm::Sp(rho()*transferRate_, eqn.psi());
        }
        else
        {
            eqn -= alpha()*rho()*transferRate_;
        }
    }
    else
    {
        FatalErrorInFunction
            << "Support for field " << rho.name() << " is not implemented"
            << exit(FatalError);
    }
}


void Foam::fv::VoFFilmTransfer::addSup
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

    if (&he == &thermo_.he())
    {
        // Explicit transfer into the VoF
        eqn +=
            filmVoFTransferRate<scalar>
            (
                &filmVoFTransfer::heTransferRate,
                dimEnergy
            );

        // Potentially implicit transfer out of the VoF
        if (&he == &eqn.psi())
        {
            eqn -= fvm::Sp(alpha()*rho()*transferRate_, eqn.psi());
        }
        else
        {
            eqn -= alpha()*rho()*he*transferRate_;
        }
    }
    else
    {
        FatalErrorInFunction
            << "Support for field " << he.name() << " is not implemented"
            << exit(FatalError);
    }
}


void Foam::fv::VoFFilmTransfer::addSup
(
    const volScalarField& rho,
    const volVectorField& U,
    fvMatrix<vector>& eqn
) const
{
    if (debug)
    {
        Info<< type() << ": applying source to " << eqn.psi().name() << endl;
    }

    if (&U == &VoF_.U)
    {
        // Explicit transfer into the VoF
        eqn +=
            filmVoFTransferRate<vector>
            (
                &filmVoFTransfer::UTransferRate,
                dimMass*dimVelocity
            );

        // Potentially implicit transfer out of the VoF
        if (&U == &eqn.psi())
        {
            eqn -= fvm::Sp(alpha_()*rho()*transferRate_, eqn.psi());
        }
        else
        {
            eqn -= alpha_()*rho()*U*transferRate_;
        }
    }
    else
    {
        FatalErrorInFunction
            << "Support for field " << U.name() << " is not implemented"
            << exit(FatalError);
    }
}


template<class Type, class FieldType>
inline Foam::tmp<Foam::Field<Type>> Foam::fv::VoFFilmTransfer::TransferRate
(
    const FieldType& f
) const
{
    const fvPatch& VoFFilmPatch = mesh().boundary()[filmPatchi_];

    const labelList& faceCells = VoFFilmPatch.faceCells();


    return tmp<Field<Type>>
    (
        new Field<Type>
        (
            UIndirectList<Type>
            (
                alpha_()*transferRate_*mesh().V()*f,
                faceCells
            )
        )
    );
}


Foam::tmp<Foam::scalarField>
Foam::fv::VoFFilmTransfer::rhoTransferRate() const
{
    return TransferRate<scalar>(thermo_.rho()());
}


Foam::tmp<Foam::scalarField>
Foam::fv::VoFFilmTransfer::heTransferRate() const
{
    return TransferRate<scalar>(thermo_.rho()()*thermo_.he()());
}


Foam::tmp<Foam::vectorField>
Foam::fv::VoFFilmTransfer::UTransferRate() const
{
    return TransferRate<vector>(thermo_.rho()()*VoF_.U());
}


void Foam::fv::VoFFilmTransfer::topoChange(const polyTopoChangeMap&)
{
    transferRate_.setSize(mesh().nCells());
}


void Foam::fv::VoFFilmTransfer::mapMesh(const polyMeshMap& map)
{
    transferRate_.setSize(mesh().nCells());
}


void Foam::fv::VoFFilmTransfer::distribute(const polyDistributionMap&)
{
    transferRate_.setSize(mesh().nCells());
}


bool Foam::fv::VoFFilmTransfer::movePoints()
{
    return true;
}


// ************************************************************************* //
