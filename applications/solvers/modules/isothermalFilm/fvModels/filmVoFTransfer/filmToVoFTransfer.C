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

#include "filmToVoFTransfer.H"
#include "VoFtoFilmTransfer.H"
#include "mappedPatchBase.H"
#include "compressibleVoF.H"
#include "fvmSup.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

namespace Foam
{
    namespace fv
    {
        defineTypeNameAndDebug(filmToVoFTransfer, 0);

        addToRunTimeSelectionTable
        (
            fvModel,
            filmToVoFTransfer,
            dictionary
        );
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::filmToVoFTransfer::filmToVoFTransfer
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

Foam::wordList Foam::fv::filmToVoFTransfer::addSupFields() const
{
    return wordList
    (
        {
            film_.alpha.name(),
            film_.thermo.he().name(),
            film_.U.name()
        }
    );
}


void Foam::fv::filmToVoFTransfer::correct()
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

    const label patchiVoF = film_.surfacePatchMap().nbrPolyPatch().index();

    const VoFtoFilmTransfer& VoFtoFilm
    (
        refCast<VoFtoFilmTransfer>(VoF_.fvModels()[0])
    );

    const scalarField alphaVoF
    (
        film_.surfacePatchMap().fromNeigbour
        (
            VoFtoFilm.alpha().boundaryField()[patchiVoF]
        )
    );

    const scalarField deltaCoeffsVoF
    (
        film_.surfacePatchMap().fromNeigbour
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
            transferRate_[celli] = -transferRateCoeff_/deltaT;
        }
    }
}


template<class Type, class TransferRateFunc>
Foam::tmp<Foam::VolInternalField<Type>>
inline Foam::fv::filmToVoFTransfer::VoFToFilmTransferRate
(
    TransferRateFunc transferRateFunc,
    const dimensionSet& dimProp
) const
{
    const Foam::fvModels& fvModels
    (
        fvModels::New
        (
            refCast<const fvMesh>(film_.surfacePatchMap().nbrMesh())
        )
    );

    const VoFtoFilmTransfer* VoFtoFilmPtr = nullptr;

    forAll(fvModels, i)
    {
        if (isType<VoFtoFilmTransfer>(fvModels[i]))
        {
            const VoFtoFilmTransfer& VoFtoFilm
            (
                refCast<const VoFtoFilmTransfer>(fvModels[i])
            );

            if
            (
                VoFtoFilm.filmPatchIndex()
             == film_.surfacePatchMap().nbrPolyPatch().index()
            )
            {
                VoFtoFilmPtr = &VoFtoFilm;
            }
        }
    }

    if (!VoFtoFilmPtr)
    {
        FatalErrorInFunction
            << "Cannot find VoFtoFilmTransfer fvModel for this film "
               "in VoF region " << film_.surfacePatchMap().nbrMesh().name()
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

    UIndirectList<Type>(tSu.ref(), film_.surfacePatch().faceCells()) =
        film_.surfacePatchMap().fromNeigbour
        (
            (VoFtoFilmPtr->*transferRateFunc)()
        );

    return tSu/mesh().V();
}


void Foam::fv::filmToVoFTransfer::addSup
(
    const volScalarField& rho,
    fvMatrix<scalar>& eqn,
    const word& fieldName
) const
{
    if (debug)
    {
        Info<< type() << ": applying source to " << eqn.psi().name() << endl;
    }

    if (fieldName == film_.alpha.name())
    {
        eqn +=
            VoFToFilmTransferRate<scalar>
            (
                &VoFtoFilmTransfer::rhoTransferRate,
                dimMass
            )
          + fvm::Sp(transferRate_*rho(), eqn.psi());
    }
    else
    {
        FatalErrorInFunction
            << "Support for field " << fieldName << " is not implemented"
            << exit(FatalError);
    }
}


void Foam::fv::filmToVoFTransfer::addSup
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

    if (fieldName == film_.thermo.he().name())
    {
        eqn +=
            VoFToFilmTransferRate<scalar>
            (
                &VoFtoFilmTransfer::heTransferRate,
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


void Foam::fv::filmToVoFTransfer::addSup
(
    const volScalarField& alpha,
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
        VoFToFilmTransferRate<vector>
        (
            &VoFtoFilmTransfer::UTransferRate,
            dimMass*dimVelocity
        )
      + fvm::Sp(alpha()*rho()*transferRate_, eqn.psi());
}


template<class Type, class FieldType>
inline Foam::tmp<Foam::Field<Type>> Foam::fv::filmToVoFTransfer::TransferRate
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
                -film_.alpha()*transferRate_*mesh().V()*f,
                faceCells
            )
        )
    );
}


Foam::tmp<Foam::scalarField>
Foam::fv::filmToVoFTransfer::transferRate() const
{
    return TransferRate<scalar>(oneField());
}


Foam::tmp<Foam::scalarField>
Foam::fv::filmToVoFTransfer::rhoTransferRate() const
{
    return TransferRate<scalar>(film_.thermo.rho()());
}


Foam::tmp<Foam::scalarField>
Foam::fv::filmToVoFTransfer::heTransferRate() const
{
    return TransferRate<scalar>(film_.thermo.rho()()*film_.thermo.he()());
}


Foam::tmp<Foam::vectorField>
Foam::fv::filmToVoFTransfer::UTransferRate() const
{
    return TransferRate<vector>(film_.thermo.rho()()*film_.U());
}


void Foam::fv::filmToVoFTransfer::topoChange(const polyTopoChangeMap&)
{
    transferRate_.setSize(mesh().nCells());
}


void Foam::fv::filmToVoFTransfer::mapMesh(const polyMeshMap& map)
{
    transferRate_.setSize(mesh().nCells());
}


void Foam::fv::filmToVoFTransfer::distribute(const polyDistributionMap&)
{
    transferRate_.setSize(mesh().nCells());
}


bool Foam::fv::filmToVoFTransfer::movePoints()
{
    return true;
}


// ************************************************************************* //
