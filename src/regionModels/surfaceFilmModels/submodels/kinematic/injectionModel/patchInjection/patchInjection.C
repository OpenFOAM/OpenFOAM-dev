/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2015-2018 OpenFOAM Foundation
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

#include "patchInjection.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace surfaceFilmModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(patchInjection, 0);
addToRunTimeSelectionTable(injectionModel, patchInjection, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

patchInjection::patchInjection
(
    surfaceFilmRegionModel& film,
    const dictionary& dict
)
:
    injectionModel(type(), film, dict),
    deltaStable_(coeffDict_.lookupOrDefault<scalar>("deltaStable", 0.0))
{
    const polyBoundaryMesh& pbm = film.regionMesh().boundaryMesh();
    patchIDs_.setSize
    (
        pbm.size() - film.regionMesh().globalData().processorPatches().size()
    );

    if (coeffDict_.found("patches"))
    {
        wordReList patchNames(coeffDict_.lookup("patches"));
        const labelHashSet patchSet = pbm.patchSet(patchNames);

        Info<< "        applying to patches:" << nl;

        label pidi = 0;
        forAllConstIter(labelHashSet, patchSet, iter)
        {
            label patchi = iter.key();
            patchIDs_[pidi++] = patchi;
            Info<< "            " << pbm[patchi].name() << endl;
        }
        patchIDs_.setSize(pidi);
        patchInjectedMasses_.setSize(pidi, 0);
    }
    else
    {
        Info<< "            applying to all patches" << endl;

        forAll(patchIDs_, patchi)
        {
            patchIDs_[patchi] = patchi;
        }

        patchInjectedMasses_.setSize(patchIDs_.size(), 0);
    }

    if (!patchIDs_.size())
    {
        FatalErrorInFunction
            << "No patches selected"
            << exit(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

patchInjection::~patchInjection()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void patchInjection::correct
(
    scalarField& availableMass,
    scalarField& massToInject,
    scalarField& diameterToInject
)
{
    // Do not correct if no patches selected
    if (!patchIDs_.size()) return;

    const scalarField& delta = film().delta();
    const scalarField& rho = film().rho();
    const scalarField& magSf = film().magSf();

    const polyBoundaryMesh& pbm = film().regionMesh().boundaryMesh();

    forAll(patchIDs_, pidi)
    {
        label patchi = patchIDs_[pidi];
        const polyPatch& pp = pbm[patchi];
        const labelList& faceCells = pp.faceCells();

        // Accumulate the total mass removed from patch
        scalar dMassPatch = 0;

        forAll(faceCells, fci)
        {
            label celli = faceCells[fci];

            scalar ddelta = max(0.0, delta[celli] - deltaStable_);
            scalar dMass = ddelta*rho[celli]*magSf[celli];
            massToInject[celli] += dMass;
            availableMass[celli] -= dMass;
            dMassPatch += dMass;
        }

        patchInjectedMasses_[pidi] += dMassPatch;
        addToInjectedMass(dMassPatch);
    }

    injectionModel::correct();

    if (writeTime())
    {
        scalarField patchInjectedMasses0
        (
            getModelProperty<scalarField>
            (
                "patchInjectedMasses",
                scalarField(patchInjectedMasses_.size(), 0)
            )
        );

        scalarField patchInjectedMassTotals(patchInjectedMasses_);
        Pstream::listCombineGather(patchInjectedMassTotals, plusEqOp<scalar>());
        patchInjectedMasses0 += patchInjectedMassTotals;

        setModelProperty<scalarField>
        (
            "patchInjectedMasses",
            patchInjectedMasses0
        );

        patchInjectedMasses_ = 0;
    }
}


void patchInjection::patchInjectedMassTotals(scalarField& patchMasses) const
{
    // Do not correct if no patches selected
    if (!patchIDs_.size()) return;

    scalarField patchInjectedMasses
    (
        getModelProperty<scalarField>
        (
            "patchInjectedMasses",
            scalarField(patchInjectedMasses_.size(), 0)
        )
    );

    scalarField patchInjectedMassTotals(patchInjectedMasses_);
    Pstream::listCombineGather(patchInjectedMassTotals, plusEqOp<scalar>());

    forAll(patchIDs_, pidi)
    {
        label patchi = patchIDs_[pidi];
        patchMasses[patchi] +=
            patchInjectedMasses[pidi] + patchInjectedMassTotals[pidi];
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace surfaceFilmModels
} // End namespace regionModels
} // End namespace Foam

// ************************************************************************* //
