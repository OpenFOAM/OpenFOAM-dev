/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2017-2022 OpenFOAM Foundation
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

#include "VoFPatchTransfer.H"
#include "compressibleTwoPhaseMixture.H"
#include "thermoSurfaceFilm.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace surfaceFilmModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(VoFPatchTransfer, 0);
addToRunTimeSelectionTable(transferModel, VoFPatchTransfer, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

VoFPatchTransfer::VoFPatchTransfer
(
    surfaceFilm& film,
    const dictionary& dict
)
:
    transferModel(type(), film, dict),
    deltaFactorToVoF_
    (
        coeffDict_.lookupOrDefault<scalar>("deltaFactorToVoF", 1.0)
    ),
    deltaFactorToFilm_
    (
        coeffDict_.lookupOrDefault<scalar>("deltaFactorToFilm", 0.5)
    ),
    alphaToVoF_
    (
        coeffDict_.lookupOrDefault<scalar>("alphaToVoF", 0.5)
    ),
    alphaToFilm_
    (
        coeffDict_.lookupOrDefault<scalar>("alphaToFilm", 0.1)
    ),
    transferRateCoeff_
    (
        coeffDict_.lookupOrDefault<scalar>("transferRateCoeff", 0.1)
    )
{
    const polyBoundaryMesh& pbm = film.mesh().boundaryMesh();
    patchIDs_.setSize
    (
        pbm.size() - film.mesh().globalData().processorPatches().size()
    );

    if (coeffDict_.found("patches"))
    {
        const wordReList patchNames(coeffDict_.lookup("patches"));
        const labelHashSet patchSet = pbm.patchSet(patchNames);

        Info<< "        applying to patches:" << nl;

        label pidi = 0;
        forAllConstIter(labelHashSet, patchSet, iter)
        {
            const label patchi = iter.key();
            patchIDs_[pidi++] = patchi;
            Info<< "            " << pbm[patchi].name() << endl;
        }
        patchIDs_.setSize(pidi);
        patchTransferredMasses_.setSize(pidi, 0);
    }
    else
    {
        Info<< "            applying to all patches" << endl;

        forAll(patchIDs_, patchi)
        {
            patchIDs_[patchi] = patchi;
        }

        patchTransferredMasses_.setSize(patchIDs_.size(), 0);
    }

    if (!patchIDs_.size())
    {
        FatalErrorInFunction
            << "No patches selected"
            << exit(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

VoFPatchTransfer::~VoFPatchTransfer()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void VoFPatchTransfer::correct
(
    scalarField& availableMass,
    scalarField& massToTransfer,
    vectorField& momentumToTransfer
)
{
    NotImplemented;
}


void VoFPatchTransfer::correct
(
    scalarField& availableMass,
    scalarField& massToTransfer,
    vectorField& momentumToTransfer,
    scalarField& energyToTransfer
)
{
    // Do not correct if no patches selected
    if (!patchIDs_.size()) return;

    const thermoSurfaceFilm& film = filmType<thermoSurfaceFilm>();

    const scalarField& delta = film.delta();
    const scalarField& rho = film.rho();
    const scalarField& magSf = film.magSf();

    const polyBoundaryMesh& pbm = film.mesh().boundaryMesh();


    const compressibleTwoPhaseMixture& thermo
    (
        film.primaryMesh().lookupObject<compressibleTwoPhaseMixture>
        (
            "phaseProperties"
        )
    );

    const volVectorField& UVoF
    (
        film.primaryMesh().lookupObject<volVectorField>("U")
    );

    const volScalarField& alphaVoF = thermo.alpha1();
    const volScalarField& rhoVoF = thermo.thermo1().rho()();
    const volScalarField& heVoF = thermo.thermo1().he();
    const volScalarField& TVoF = thermo.thermo1().T();
    const volScalarField CpVoF(thermo.thermo1().Cp());

    forAll(patchIDs_, pidi)
    {
        const label patchi = patchIDs_[pidi];
        label primaryPatchi = -1;

        forAll(film.intCoupledPatchIDs(), i)
        {
            const label filmPatchi = film.intCoupledPatchIDs()[i];

            if (filmPatchi == patchi)
            {
                primaryPatchi = film.primaryPatchIDs()[i];
            }
        }

        if (primaryPatchi != -1)
        {
            scalarField deltaCoeffs
            (
                film.primaryMesh().boundary()[primaryPatchi].deltaCoeffs()
            );
            film.toRegion(patchi, deltaCoeffs);

            scalarField alphap(alphaVoF.boundaryField()[primaryPatchi]);
            film.toRegion(patchi, alphap);

            scalarField rhop(rhoVoF.boundaryField()[primaryPatchi]);
            film.toRegion(patchi, rhop);

            vectorField Up(UVoF.boundaryField()[primaryPatchi]);
            film.toRegion(patchi, Up);

            scalarField hp(heVoF.boundaryField()[primaryPatchi]);
            film.toRegion(patchi, hp);

            scalarField Tp(TVoF.boundaryField()[primaryPatchi]);
            film.toRegion(patchi, Tp);

            scalarField Cpp(CpVoF.boundaryField()[primaryPatchi]);
            film.toRegion(patchi, Cpp);

            scalarField Vp
            (
                film.primaryMesh().boundary()[primaryPatchi]
               .patchInternalField(film.primaryMesh().V())
            );
            film.toRegion(patchi, Vp);

            const polyPatch& pp = pbm[patchi];
            const labelList& faceCells = pp.faceCells();

            const vectorField& U = film.U();
            const scalarField& he = film.thermo().he();

            // Accumulate the total mass removed from patch
            scalar dMassPatch = 0;

            forAll(faceCells, facei)
            {
                const label celli = faceCells[facei];

                scalar dMass = 0;

                if
                (
                    delta[celli] > 2*deltaFactorToVoF_/deltaCoeffs[facei]
                 || alphap[facei] > alphaToVoF_
                )
                {
                    dMass =
                        transferRateCoeff_*delta[celli]*rho[celli]*magSf[celli];

                    massToTransfer[celli] += dMass;
                    momentumToTransfer[celli] += dMass*U[celli];
                    energyToTransfer[celli] += dMass*he[celli];
                }

                if
                (
                    alphap[facei] > 0
                 && delta[celli] < 2*deltaFactorToFilm_/deltaCoeffs[facei]
                 && alphap[facei] < alphaToFilm_
                )
                {
                    dMass =
                        -transferRateCoeff_*alphap[facei]*rhop[facei]*Vp[facei];

                    massToTransfer[celli] += dMass;
                    momentumToTransfer[celli] += dMass*Up[facei];
                    energyToTransfer[celli] += dMass*hp[facei];
                }

                availableMass[celli] -= dMass;
                dMassPatch += dMass;
            }

            patchTransferredMasses_[pidi] += dMassPatch;
            addToTransferredMass(dMassPatch);
        }
    }

    transferModel::correct();

    if (writeTime())
    {
        scalarField patchTransferredMasses0
        (
            getModelProperty<scalarField>
            (
                "patchTransferredMasses",
                scalarField(patchTransferredMasses_.size(), 0)
            )
        );

        scalarField patchTransferredMassTotals(patchTransferredMasses_);
        Pstream::listCombineGather
        (
            patchTransferredMassTotals,
            plusEqOp<scalar>()
        );
        patchTransferredMasses0 += patchTransferredMassTotals;

        setModelProperty<scalarField>
        (
            "patchTransferredMasses",
            patchTransferredMasses0
        );

        patchTransferredMasses_ = 0;
    }
}


void VoFPatchTransfer::patchTransferredMassTotals
(
    scalarField& patchMasses
) const
{
    // Do not correct if no patches selected
    if (!patchIDs_.size()) return;

    scalarField patchTransferredMasses
    (
        getModelProperty<scalarField>
        (
            "patchTransferredMasses",
            scalarField(patchTransferredMasses_.size(), 0)
        )
    );

    scalarField patchTransferredMassTotals(patchTransferredMasses_);
    Pstream::listCombineGather(patchTransferredMassTotals, plusEqOp<scalar>());

    forAll(patchIDs_, pidi)
    {
        const label patchi = patchIDs_[pidi];
        patchMasses[patchi] +=
            patchTransferredMasses[pidi] + patchTransferredMassTotals[pidi];
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace surfaceFilmModels
} // End namespace Foam

// ************************************************************************* //
