/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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

#include "removeInjection.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace surfaceFilmModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(removeInjection, 0);
addToRunTimeSelectionTable(injectionModel, removeInjection, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

removeInjection::removeInjection
(
    surfaceFilmModel& owner,
    const dictionary& dict
)
:
    injectionModel(type(), owner, dict),
    deltaStable_(coeffDict_.lookupOrDefault<scalar>("deltaStable", 0.0)),
    mask_(owner.regionMesh().nCells(), -1)
{
    wordReList patches;
    if (coeffDict_.readIfPresent("patches", patches))
    {
        Info<< "        applying to patches:" << nl;
        const polyBoundaryMesh& pbm = owner.regionMesh().boundaryMesh();
        const labelHashSet patchSet = pbm.patchSet(patches);

        forAllConstIter(labelHashSet, patchSet, iter)
        {
            label patchI = iter.key();
            const polyPatch& pp = pbm[patchI];
            UIndirectList<scalar>(mask_, pp.faceCells()) = 1.0;

            Info<< "            " << pp.name() << endl;
        }
    }
    else
    {
        Info<< "            applying to all patches" << endl;
        mask_ = 1.0;
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

removeInjection::~removeInjection()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void removeInjection::correct
(
    scalarField& availableMass,
    scalarField& massToInject,
    scalarField& diameterToInject
)
{
    const scalarField& delta = owner().delta();
    const scalarField& rho = owner().rho();
    const scalarField& magSf = owner().magSf();

    forAll(delta, cellI)
    {
        if (mask_[cellI] > 0)
        {
            scalar ddelta = max(0.0, delta[cellI] - deltaStable_);
            scalar dMass = ddelta*rho[cellI]*magSf[cellI];
            massToInject[cellI] += dMass;
            availableMass[cellI] -= dMass;

            addToInjectedMass(dMass);
        }
    }

    injectionModel::correct();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace surfaceFilmModels
} // End namespace regionModels
} // End namespace Foam

// ************************************************************************* //
