/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2017-2019 OpenFOAM Foundation
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
#include "twoPhaseMixtureThermo.H"
#include "thermoSingleLayer.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace surfaceFilmModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(VoFPatchTransfer, 0);
addToRunTimeSelectionTable(transferModel, VoFPatchTransfer, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

VoFPatchTransfer::VoFPatchTransfer
(
    surfaceFilmRegionModel& film,
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
    const polyBoundaryMesh& pbm = film.regionMesh().boundaryMesh();
    patchIDs_.setSize
    (
        pbm.size() - film.regionMesh().globalData().processorPatches().size()
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
    scalarField& massToTransfer
)
{
    NotImplemented;
}


void VoFPatchTransfer::correct
(
    scalarField& availableMass,
    scalarField& massToTransfer,
    scalarField& energyToTransfer
)
{
    // Do not correct if no patches selected
    if (!patchIDs_.size()) return;

    const thermoSingleLayer& film = filmType<thermoSingleLayer>();

    const scalarField& delta = film.delta();
    const scalarField& rho = film.rho();
    const scalarField& magSf = film.magSf();

    const polyBoundaryMesh& pbm = film.regionMesh().boundaryMesh();


    const twoPhaseMixtureThermo& thermo
    (
        film.primaryMesh().lookupObject<twoPhaseMixtureThermo>
        (
            twoPhaseMixtureThermo::dictName
        )
    );

    const volScalarField& heVoF = thermo.thermo1().he();
    const volScalarField& TVoF = thermo.thermo1().T();
    const volScalarField CpVoF(thermo.thermo1().Cp());
    const volScalarField& rhoVoF = thermo.thermo1().rho()();
    const volScalarField& alphaVoF = thermo.alpha1();

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

            scalarField hp(heVoF.boundaryField()[primaryPatchi]);
            film.toRegion(patchi, hp);

            scalarField Tp(TVoF.boundaryField()[primaryPatchi]);
            film.toRegion(patchi, Tp);

            scalarField Cpp(CpVoF.boundaryField()[primaryPatchi]);
            film.toRegion(patchi, Cpp);

            scalarField rhop(rhoVoF.boundaryField()[primaryPatchi]);
            film.toRegion(patchi, rhop);

            scalarField alphap(alphaVoF.boundaryField()[primaryPatchi]);
            film.toRegion(patchi, alphap);

            scalarField Vp
            (
                film.primaryMesh().boundary()[primaryPatchi]
               .patchInternalField(film.primaryMesh().V())
            );
            film.toRegion(patchi, Vp);

            const polyPatch& pp = pbm[patchi];
            const labelList& faceCells = pp.faceCells();

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
                    energyToTransfer[celli] += dMass*film.h()[celli];
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
} // End namespace regionModels
} // End namespace Foam

// ************************************************************************* //
