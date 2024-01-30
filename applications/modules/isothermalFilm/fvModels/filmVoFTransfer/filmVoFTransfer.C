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

#include "filmVoFTransfer.H"
#include "VoFFilmTransfer.H"
#include "mappedFvPatchBaseBase.H"
#include "compressibleVoF.H"
#include "fvmSup.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

namespace Foam
{
    namespace fv
    {
        defineTypeNameAndDebug(filmVoFTransfer, 0);

        addToRunTimeSelectionTable
        (
            fvModel,
            filmVoFTransfer,
            dictionary
        );
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::filmVoFTransfer::filmVoFTransfer
(
    const word& sourceName,
    const word& modelType,
    const fvMesh& mesh,
    const dictionary& dict
)
:
    fvModel(sourceName, modelType, mesh, dict),
    film_(mesh.lookupObject<solvers::isothermalFilm>(solver::typeName)),
    curTimeIndex_(-1),
    deltaFactorToVoF_
    (
        dict.lookupOrDefault<scalar>("deltaFactorToVoF", 1.0)
    ),
    alphaToVoF_
    (
        dict.lookupOrDefault<scalar>("alphaToVoF", 0.5)
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

Foam::wordList Foam::fv::filmVoFTransfer::addSupFields() const
{
    return wordList
    {
        film_.alpha.name(),
        film_.thermo.he().name(),
        film_.U.name()
    };
}


void Foam::fv::filmVoFTransfer::correct()
{
    if (curTimeIndex_ == mesh().time().timeIndex())
    {
        return;
    }

    curTimeIndex_ = mesh().time().timeIndex();

    const scalar deltaT = mesh().time().deltaTValue();


    // Film properties

    const labelList& faceCells = film_.surfacePatch().faceCells();
    const scalarField& delta = film_.delta();


    // VoF properties

    const solvers::compressibleVoF& VoF_
    (
        film_.surfacePatchMap().nbrMesh().lookupObject<solvers::compressibleVoF>
        (
            solver::typeName
        )
    );

    const label patchiVoF = film_.surfacePatchMap().nbrFvPatch().index();

    const VoFFilmTransfer& VoFFilm(this->VoFFilm(VoF_.fvModels()));

    const scalarField alphaVoF
    (
        film_.surfacePatchMap().fromNeighbour
        (
            VoFFilm.alpha().boundaryField()[patchiVoF]
        )
    );

    const scalarField deltaCoeffsVoF
    (
        film_.surfacePatchMap().fromNeighbour
        (
            VoF_.mesh.boundary()[patchiVoF].deltaCoeffs()
        )
    );

    // Reset the transfer rate
    transferRate_ = Zero;

    forAll(faceCells, facei)
    {
        const label celli = faceCells[facei];

        if
        (
            delta[celli] > 2*deltaFactorToVoF_/deltaCoeffsVoF[facei]
         || alphaVoF[facei] > alphaToVoF_
        )
        {
            transferRate_[celli] = transferRateCoeff_/deltaT;
        }
    }
}


const Foam::fv::VoFFilmTransfer& Foam::fv::filmVoFTransfer::VoFFilm
(
    const Foam::fvModels& fvModels
) const
{
    const VoFFilmTransfer* VoFFilmPtr = nullptr;

    forAll(fvModels, i)
    {
        if (isType<VoFFilmTransfer>(fvModels[i]))
        {
            const VoFFilmTransfer& VoFFilm
            (
                refCast<const VoFFilmTransfer>(fvModels[i])
            );

            if
            (
                VoFFilm.filmPatchIndex()
             == film_.surfacePatchMap().nbrFvPatch().index()
            )
            {
                VoFFilmPtr = &VoFFilm;
            }
        }
    }

    if (!VoFFilmPtr)
    {
        FatalErrorInFunction
            << "Cannot find VoFFilmTransfer fvModel for this film "
               "in VoF region " << film_.surfacePatchMap().nbrMesh().name()
            << exit(FatalError);
    }

    return *VoFFilmPtr;
}


template<class Type, class TransferRateFunc>
Foam::tmp<Foam::VolInternalField<Type>>
inline Foam::fv::filmVoFTransfer::VoFToFilmTransferRate
(
    TransferRateFunc transferRateFunc,
    const dimensionSet& dimProp
) const
{
    const Foam::fvModels& fvModels =
        fvModels::New(film_.surfacePatchMap().nbrMesh());

    tmp<VolInternalField<Type>> tSu
    (
        VolInternalField<Type>::New
        (
            "Su",
            mesh(),
            dimensioned<Type>(dimProp/dimTime, Zero)
        )
    );

    UIndirectList<Type>(tSu.ref(), film_.surfacePatch().faceCells()) =
        film_.surfacePatchMap().fromNeighbour
        (
            (VoFFilm(fvModels).*transferRateFunc)()
        );

    return tSu/mesh().V();
}


void Foam::fv::filmVoFTransfer::addSup
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

    if (&alpha == &film_.alpha)
    {
        // Explicit transfer into the film
        eqn +=
            VoFToFilmTransferRate<scalar>
            (
                &VoFFilmTransfer::rhoTransferRate,
                dimMass
            );

        // Potentially implicit transfer out of the film
        if (&alpha == &eqn.psi())
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
            << "Support for field " << alpha.name() << " is not implemented"
            << exit(FatalError);
    }
}


void Foam::fv::filmVoFTransfer::addSup
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

    if (&he == &film_.thermo.he())
    {
        // Explicit transfer into the film
        eqn +=
            VoFToFilmTransferRate<scalar>
            (
                &VoFFilmTransfer::heTransferRate,
                dimEnergy
            );

        // Potentially implicit transfer out of the film
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


void Foam::fv::filmVoFTransfer::addSup
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

    if (&U == &film_.U)
    {
        // Explicit transfer into the film
        eqn +=
            VoFToFilmTransferRate<vector>
            (
                &VoFFilmTransfer::UTransferRate,
                dimMomentum
            );

        // Potentially implicit transfer out of the film
        if (&U == &eqn.psi())
        {
            eqn -= fvm::Sp(alpha()*rho()*transferRate_, eqn.psi());
        }
        else
        {
            eqn -= alpha()*rho()*U()*transferRate_;
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
inline Foam::tmp<Foam::Field<Type>> Foam::fv::filmVoFTransfer::TransferRate
(
    const FieldType& f
) const
{
    const labelList& faceCells = film_.surfacePatch().faceCells();

    return tmp<Field<Type>>
    (
        new Field<Type>
        (
            UIndirectList<Type>
            (
                film_.alpha()*transferRate_*mesh().V()*f,
                faceCells
            )
        )
    );
}


Foam::tmp<Foam::scalarField>
Foam::fv::filmVoFTransfer::transferRate() const
{
    return TransferRate<scalar>(oneField());
}


Foam::tmp<Foam::scalarField>
Foam::fv::filmVoFTransfer::rhoTransferRate() const
{
    return TransferRate<scalar>(film_.thermo.rho()());
}


Foam::tmp<Foam::scalarField>
Foam::fv::filmVoFTransfer::heTransferRate() const
{
    return TransferRate<scalar>(film_.thermo.rho()()*film_.thermo.he()());
}


Foam::tmp<Foam::vectorField>
Foam::fv::filmVoFTransfer::UTransferRate() const
{
    return TransferRate<vector>(film_.thermo.rho()()*film_.U());
}


void Foam::fv::filmVoFTransfer::topoChange(const polyTopoChangeMap&)
{
    transferRate_.setSize(mesh().nCells());
}


void Foam::fv::filmVoFTransfer::mapMesh(const polyMeshMap& map)
{
    transferRate_.setSize(mesh().nCells());
}


void Foam::fv::filmVoFTransfer::distribute(const polyDistributionMap&)
{
    transferRate_.setSize(mesh().nCells());
}


bool Foam::fv::filmVoFTransfer::movePoints()
{
    return true;
}


// ************************************************************************* //
