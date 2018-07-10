/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2017-2018 OpenFOAM Foundation
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

#include "transferModelList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace surfaceFilmModels
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

transferModelList::transferModelList(surfaceFilmRegionModel& film)
:
    PtrList<transferModel>(),
    filmSubModelBase(film)
{}


transferModelList::transferModelList
(
    surfaceFilmRegionModel& film,
    const dictionary& dict
)
:
    PtrList<transferModel>(),
    filmSubModelBase
    (
        "transferModelList",
        film,
        dict,
        "transferModelList",
        "transferModelList"
    ),
    massTransferred_(film.intCoupledPatchIDs().size(), 0.0)
{
    const wordList activeModels
    (
        dict.lookupOrDefault("transferModels", wordList())
    );

    wordHashSet models;
    forAll(activeModels, i)
    {
        models.insert(activeModels[i]);
    }

    Info<< "    Selecting film transfer models" << endl;
    if (models.size() > 0)
    {
        this->setSize(models.size());

        label i = 0;
        forAllConstIter(wordHashSet, models, iter)
        {
            const word& model = iter.key();
            set(i, transferModel::New(film, dict, model));
            i++;
        }
    }
    else
    {
        Info<< "        none" << endl;
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

transferModelList::~transferModelList()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void transferModelList::correct
(
    scalarField& availableMass,
    volScalarField& massToTransfer
)
{
    // Correct models that accumulate mass and diameter transfers
    forAll(*this, i)
    {
        operator[](i).correct(availableMass, massToTransfer);
    }

    // Push values to boundaries ready for transfer to the primary region
    massToTransfer.correctBoundaryConditions();

    const labelList& patchIDs = film().intCoupledPatchIDs();

    forAll(patchIDs, i)
    {
        label patchi = patchIDs[i];
        massTransferred_[i] =
            massTransferred_[i] + sum(massToTransfer.boundaryField()[patchi]);
    }
}


void transferModelList::correct
(
    scalarField& availableMass,
    volScalarField& massToTransfer,
    volScalarField& energyToTransfer
)
{
    // Correct models that accumulate mass and diameter transfers
    forAll(*this, i)
    {
        operator[](i).correct(availableMass, massToTransfer, energyToTransfer);
    }

    // Push values to boundaries ready for transfer to the primary region
    massToTransfer.correctBoundaryConditions();
    energyToTransfer.correctBoundaryConditions();

    const labelList& patchIDs = film().intCoupledPatchIDs();

    forAll(patchIDs, i)
    {
        label patchi = patchIDs[i];
        massTransferred_[i] =
            massTransferred_[i] + sum(massToTransfer.boundaryField()[patchi]);
    }
}


void transferModelList::info(Ostream& os)
{
    const polyBoundaryMesh& pbm = film().regionMesh().boundaryMesh();

    scalar transferredMass = 0;
    scalarField patchTransferredMasses
    (
        pbm.size() - film().regionMesh().globalData().processorPatches().size(),
        0
    );

    forAll(*this, i)
    {
        const transferModel& im = operator[](i);
        transferredMass += im.transferredMassTotal();
        im.patchTransferredMassTotals(patchTransferredMasses);
    }

    os  << indent << "transferred mass      = " << transferredMass << nl;

    forAll(patchTransferredMasses, patchi)
    {
        if (mag(patchTransferredMasses[patchi]) > vSmall)
        {
            os  << indent << indent << "from patch " << pbm[patchi].name()
                << " = " << patchTransferredMasses[patchi] << nl;
        }
    }

    scalarField mass0(massTransferred_.size(), 0);
    this->getBaseProperty("massTransferred", mass0);

    scalarField mass(massTransferred_);
    Pstream::listCombineGather(mass, plusEqOp<scalar>());
    mass += mass0;

    const labelList& patchIDs = film().intCoupledPatchIDs();

    forAll(patchIDs, i)
    {
        label patchi = patchIDs[i];
        Info<< indent << "  - patch: " << pbm[patchi].name() << ": "
            << mass[i] << endl;
    }

    if (film().time().writeTime())
    {
        setBaseProperty("massTransferred", mass);
        massTransferred_ = 0.0;
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace surfaceFilmModels
} // End namespace regionModels
} // End namespace Foam

// ************************************************************************* //
