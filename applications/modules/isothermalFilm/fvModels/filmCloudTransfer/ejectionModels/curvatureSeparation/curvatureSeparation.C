/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2023 OpenFOAM Foundation
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

#include "curvatureSeparation.H"
#include "fvcGrad.H"
#include "zeroGradientFvPatchFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace filmEjectionModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(curvatureSeparation, 0);
addToRunTimeSelectionTable
(
    ejectionModel,
    curvatureSeparation,
    dictionary
);

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

tmp<scalarField> curvatureSeparation::calcInvR1
(
    const volVectorField& U
) const
{
    const vectorField UHat(U.field()/(mag(U.field()) + rootVSmall));

    tmp<scalarField> tinvR1(-(UHat & (UHat & gradNHat_)));
    scalarField& invR1 = tinvR1.ref();

    const scalar rMin = 1e-6;
    const fvMesh& mesh = film_.mesh;
    const polyBoundaryMesh& pbm = mesh.boundaryMesh();

    forAll(patchRadii_, i)
    {
        const label patchi = patchRadii_[i].first();
        const scalar definedInvR1 = 1/max(rMin, patchRadii_[i].second());
        UIndirectList<scalar>(invR1, pbm[patchi].faceCells()) = definedInvR1;
    }

    // Filter out large radii
    const scalar rMax = 1e6;

    forAll(invR1, i)
    {
        if (mag(invR1[i]) < 1/rMax)
        {
            invR1[i] = -1;
        }
    }

    return tinvR1;
}


tmp<scalarField> curvatureSeparation::calcCosAngle
(
    const surfaceScalarField& alphaRhoPhi
) const
{
    const fvMesh& mesh = film_.mesh;

    const scalar magg(mag(film_.g.value()));
    const vector gHat(film_.g.value()/magg);

    const vectorField nf(mesh.Sf()/mesh.magSf());
    const labelUList& own = mesh.owner();
    const labelUList& nbr = mesh.neighbour();

    scalarField alphaRhoPhiMax(mesh.nCells(), -great);
    tmp<scalarField> tcosAngle
    (
        new scalarField(mesh.nCells(), 0)
    );
    scalarField& cosAngle = tcosAngle.ref();

    forAll(nbr, facei)
    {
        const label cellO = own[facei];
        const label cellN = nbr[facei];

        if (alphaRhoPhi[facei] > alphaRhoPhiMax[cellO])
        {
            alphaRhoPhiMax[cellO] = alphaRhoPhi[facei];
            cosAngle[cellO] = -gHat & nf[facei];
        }
        if (-alphaRhoPhi[facei] > alphaRhoPhiMax[cellN])
        {
            alphaRhoPhiMax[cellN] = -alphaRhoPhi[facei];
            cosAngle[cellN] = -gHat & -nf[facei];
        }
    }

    forAll(alphaRhoPhi.boundaryField(), patchi)
    {
        const fvsPatchScalarField& alphaRhoPhip =
            alphaRhoPhi.boundaryField()[patchi];

        const fvPatch& pp = alphaRhoPhip.patch();

        if (pp.coupled())
        {
            const labelList& faceCells = pp.faceCells();
            const vectorField nf(pp.nf());

            forAll(alphaRhoPhip, i)
            {
                const label celli = faceCells[i];

                if (alphaRhoPhip[i] > alphaRhoPhiMax[celli])
                {
                    alphaRhoPhiMax[celli] = alphaRhoPhip[i];
                    cosAngle[celli] = -gHat & nf[i];
                }
            }
        }
    }

    return tcosAngle;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

curvatureSeparation::curvatureSeparation
(
    const dictionary& dict,
    const solvers::isothermalFilm& film
)
:
    ejectionModel(dict, film),
    gradNHat_(fvc::grad(film_.nHat)),
    deltaByR1Min_
    (
        dict.optionalSubDict(typeName + "Coeffs")
       .lookupOrDefault<scalar>("deltaByR1Min", scalar(0))
    ),
    deltaStable_
    (
        dict.optionalSubDict(typeName + "Coeffs")
      .lookupOrDefault("deltaStable", scalar(0))
    )
{
    const List<Tuple2<word, scalar>> prIn
    (
        dict.lookupOrDefault("patchRadii", List<Tuple2<word, scalar>>::null())
    );
    const wordList& allPatchNames = film_.mesh.boundaryMesh().names();

    DynamicList<Tuple2<label, scalar>> prData(allPatchNames.size());

    labelHashSet uniquePatchIDs;

    forAllReverse(prIn, i)
    {
        labelList patchIDs = findStrings(prIn[i].first(), allPatchNames);
        forAll(patchIDs, j)
        {
            const label patchi = patchIDs[j];

            if (!uniquePatchIDs.found(patchi))
            {
                const scalar radius = prIn[i].second();
                prData.append(Tuple2<label, scalar>(patchi, radius));

                uniquePatchIDs.insert(patchi);
            }
        }
    }

    patchRadii_.transfer(prData);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

curvatureSeparation::~curvatureSeparation()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void curvatureSeparation::correct()
{
    const scalarField& delta = film_.delta();
    const scalarField& rho = film_.rho();
    const vectorField& U = film_.U();

    const tmp<volScalarField> tsigma = film_.sigma();
    const volScalarField::Internal& sigma = tsigma();

    const scalarField invR1(calcInvR1(film_.U));
    const scalarField cosAngle(calcCosAngle(film_.alphaRhoPhi));

    const scalar magg(mag(film_.g.value()));

    const scalar deltaT = film_.mesh.time().deltaTValue();

    // Calculate force balance
    const scalar Fthreshold = 1e-10;

    forAll(delta, celli)
    {
        rate_[celli] = 0;
        diameter_[celli] = 0;

        if
        (
            delta[celli] > deltaStable_
         && invR1[celli] > 0
         && delta[celli]*invR1[celli] > deltaByR1Min_
        )
        {
            const scalar R1 = 1/(invR1[celli] + rootVSmall);
            const scalar R2 = R1 + delta[celli];

            // Inertial force
            const scalar Fi =
                -delta[celli]*rho[celli]
                *magSqr(U[celli])*72/60*invR1[celli];

            // Body force
            const scalar Fb =
                -0.5*rho[celli]*magg
                *invR1[celli]*(sqr(R1) - sqr(R2))*cosAngle[celli];

            // Surface tension force
            const scalar Fs = sigma[celli]/R2;

            if (Fi + Fb + Fs + Fthreshold < 0)
            {
                // Calculate rate of separation
                rate_[celli] =
                    (delta[celli] - deltaStable_)/(deltaT*delta[celli]);

                // Set the separated droplet diameter to the film thickness
                diameter_[celli] = delta[celli];
            }
        }
    }
}


void curvatureSeparation::topoChange(const polyTopoChangeMap& map)
{
    ejectionModel::topoChange(map);
    gradNHat_ = fvc::grad(film_.nHat);
}


void curvatureSeparation::mapMesh(const polyMeshMap& map)
{
    ejectionModel::mapMesh(map);
    gradNHat_ = fvc::grad(film_.nHat);
}


void curvatureSeparation::distribute(const polyDistributionMap& map)
{
    ejectionModel::distribute(map);
    gradNHat_ = fvc::grad(film_.nHat);
}


bool curvatureSeparation::movePoints()
{
    ejectionModel::movePoints();
    gradNHat_ = fvc::grad(film_.nHat);

    return true;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace filmEjectionModels
} // End namespace Foam

// ************************************************************************* //
