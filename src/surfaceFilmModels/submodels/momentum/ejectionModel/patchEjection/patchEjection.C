/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2015-2022 OpenFOAM Foundation
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

#include "patchEjection.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace surfaceFilmSubModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(patchEjection, 0);
addToRunTimeSelectionTable(ejectionModel, patchEjection, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

patchEjection::patchEjection
(
    surfaceFilm& film,
    const dictionary& dict
)
:
    ejectionModel(type(), film, dict),
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
        patchEjectedMasses_.setSize(pidi, 0);
    }
    else
    {
        Info<< "            applying to all patches" << endl;

        forAll(patchIDs_, patchi)
        {
            patchIDs_[patchi] = patchi;
        }

        patchEjectedMasses_.setSize(patchIDs_.size(), 0);
    }

    if (!patchIDs_.size())
    {
        FatalErrorInFunction
            << "No patches selected"
            << exit(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

patchEjection::~patchEjection()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void patchEjection::correct
(
    scalarField& availableMass,
    scalarField& massToEject,
    scalarField& diameterToEject
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
            massToEject[celli] += dMass;
            availableMass[celli] -= dMass;
            dMassPatch += dMass;
        }

        patchEjectedMasses_[pidi] += dMassPatch;
        addToEjectedMass(dMassPatch);
    }

    ejectionModel::correct();

    if (writeTime())
    {
        scalarField patchEjectedMasses0
        (
            getModelProperty<scalarField>
            (
                "patchEjectedMasses",
                scalarField(patchEjectedMasses_.size(), 0)
            )
        );

        scalarField patchEjectedMassTotals(patchEjectedMasses_);
        Pstream::listCombineGather(patchEjectedMassTotals, plusEqOp<scalar>());
        patchEjectedMasses0 += patchEjectedMassTotals;

        setModelProperty<scalarField>
        (
            "patchEjectedMasses",
            patchEjectedMasses0
        );

        patchEjectedMasses_ = 0;
    }
}


void patchEjection::patchEjectedMassTotals(scalarField& patchMasses) const
{
    // Do not correct if no patches selected
    if (!patchIDs_.size()) return;

    scalarField patchEjectedMasses
    (
        getModelProperty<scalarField>
        (
            "patchEjectedMasses",
            scalarField(patchEjectedMasses_.size(), 0)
        )
    );

    scalarField patchEjectedMassTotals(patchEjectedMasses_);
    Pstream::listCombineGather(patchEjectedMassTotals, plusEqOp<scalar>());

    forAll(patchIDs_, pidi)
    {
        label patchi = patchIDs_[pidi];
        patchMasses[patchi] +=
            patchEjectedMasses[pidi] + patchEjectedMassTotals[pidi];
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace surfaceFilmSubModels
} // End namespace regionModels
} // End namespace Foam

// ************************************************************************* //
