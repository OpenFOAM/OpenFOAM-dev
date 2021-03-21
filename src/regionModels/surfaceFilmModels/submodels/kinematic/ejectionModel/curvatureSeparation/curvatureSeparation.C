/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2021 OpenFOAM Foundation
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
#include "kinematicSingleLayer.H"
#include "fvcGrad.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace surfaceFilmModels
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

tmp<volScalarField> curvatureSeparation::calcInvR1
(
    const volVectorField& U
) const
{
    const dimensionedScalar smallU("smallU", dimVelocity, rootVSmall);
    const volVectorField UHat(U/(mag(U) + smallU));

    tmp<volScalarField> tinvR1
    (
        volScalarField::New("invR1", UHat & (UHat & gradNHat_))
    );

    scalarField& invR1 = tinvR1.ref().primitiveFieldRef();

    // Apply defined patch radii
    const scalar rMin = 1e-6;
    const fvMesh& mesh = film().regionMesh();
    const polyBoundaryMesh& pbm = mesh.boundaryMesh();

    forAll(definedPatchRadii_, i)
    {
        const label patchi = definedPatchRadii_[i].first();
        const scalar definedInvR1 =
            1.0/max(rMin, definedPatchRadii_[i].second());
        UIndirectList<scalar>(invR1, pbm[patchi].faceCells()) = definedInvR1;
    }

    // Filter out large radii
    const scalar rMax = 1e6;
    forAll(invR1, i)
    {
        if (mag(invR1[i]) < 1/rMax)
        {
            invR1[i] = -1.0;
        }
    }

    if (debug && mesh.time().writeTime())
    {
        tinvR1().write();
    }

    return tinvR1;
}


tmp<scalarField> curvatureSeparation::calcCosAngle
(
    const surfaceScalarField& phi
) const
{
    const fvMesh& mesh = film().regionMesh();
    const vectorField nf(mesh.Sf()/mesh.magSf());
    const unallocLabelList& own = mesh.owner();
    const unallocLabelList& nbr = mesh.neighbour();

    scalarField phiMax(mesh.nCells(), -great);
    scalarField cosAngle(mesh.nCells(), 0);

    forAll(nbr, facei)
    {
        const label cellO = own[facei];
        const label cellN = nbr[facei];

        if (phi[facei] > phiMax[cellO])
        {
            phiMax[cellO] = phi[facei];
            cosAngle[cellO] = -gHat_ & nf[facei];
        }
        if (-phi[facei] > phiMax[cellN])
        {
            phiMax[cellN] = -phi[facei];
            cosAngle[cellN] = -gHat_ & -nf[facei];
        }
    }

    forAll(phi.boundaryField(), patchi)
    {
        const fvsPatchScalarField& phip = phi.boundaryField()[patchi];
        const fvPatch& pp = phip.patch();
        const labelList& faceCells = pp.faceCells();
        const vectorField nf(pp.nf());
        forAll(phip, i)
        {
            const label celli = faceCells[i];
            if (phip[i] > phiMax[celli])
            {
                phiMax[celli] = phip[i];
                cosAngle[celli] = -gHat_ & nf[i];
            }
        }
    }

    if (debug && mesh.time().writeTime())
    {
        volScalarField volCosAngle
        (
            IOobject
            (
                "cosAngle",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ
            ),
            mesh,
            dimensionedScalar(dimless, 0),
            zeroGradientFvPatchScalarField::typeName
        );
        volCosAngle.primitiveFieldRef() = cosAngle;
        volCosAngle.correctBoundaryConditions();
        volCosAngle.write();
    }

    return max(min(cosAngle, scalar(1)), scalar(-1));
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

curvatureSeparation::curvatureSeparation
(
    surfaceFilmRegionModel& film,
    const dictionary& dict
)
:
    ejectionModel(type(), film, dict),
    gradNHat_(fvc::grad(film.nHat())),
    deltaByR1Min_
    (
        coeffDict_.lookupOrDefault<scalar>("deltaByR1Min", scalar(0))
    ),
    deltaStable_(coeffDict_.lookupOrDefault("deltaStable", scalar(0))),
    definedPatchRadii_(),
    magG_(mag(film.g().value())),
    gHat_(Zero)
{
    if (magG_ < rootVSmall)
    {
        FatalErrorInFunction
            << "Acceleration due to gravity must be non-zero"
            << exit(FatalError);
    }

    gHat_ = film.g().value()/magG_;

    const List<Tuple2<word, scalar>> prIn
    (
        coeffDict_.lookup("definedPatchRadii")
    );
    const wordList& allPatchNames = film.regionMesh().boundaryMesh().names();

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

    definedPatchRadii_.transfer(prData);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

curvatureSeparation::~curvatureSeparation()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void curvatureSeparation::correct
(
    scalarField& availableMass,
    scalarField& massToEject,
    scalarField& diameterToEject
)
{
    const kinematicSingleLayer& film =
        refCast<const kinematicSingleLayer>(this->film());
    const fvMesh& mesh = film.regionMesh();

    const volScalarField& delta = film.delta();
    const volVectorField& U = film.U();
    const surfaceScalarField& phi = film.phi();
    const volScalarField& rho = film.rho();
    const scalarField magSqrU(magSqr(film.U()));

    const tmp<volScalarField> tsigma = film.sigma();
    const volScalarField::Internal& sigma = tsigma();

    const scalarField invR1(calcInvR1(U));
    const scalarField cosAngle(calcCosAngle(phi));

    // Calculate force balance
    const scalar Fthreshold = 1e-10;
    scalarField Fnet(mesh.nCells(), 0);
    scalarField separated(mesh.nCells(), 0);

    forAll(invR1, i)
    {
        if ((invR1[i] > 0) && (delta[i]*invR1[i] > deltaByR1Min_))
        {
            const scalar R1 = 1.0/(invR1[i] + rootVSmall);
            const scalar R2 = R1 + delta[i];

            // Inertial force
            const scalar Fi = -delta[i]*rho[i]*magSqrU[i]*72.0/60.0*invR1[i];

            // Body force
            const scalar Fb =
              - 0.5*rho[i]*magG_*invR1[i]*(sqr(R1) - sqr(R2))*cosAngle[i];

            // Surface force
            const scalar Fs = sigma[i]/R2;

            Fnet[i] = Fi + Fb + Fs;

            if (Fnet[i] + Fthreshold < 0 && delta[i] > deltaStable_)
            {
                separated[i] = 1.0;
            }
        }
    }

    // Eject all available mass
    massToEject = separated*availableMass;
    diameterToEject = separated*delta;
    availableMass -= separated*availableMass;

    addToEjectedMass(sum(separated*availableMass));

    if (debug && mesh.time().writeTime())
    {
        volScalarField volFnet
        (
            IOobject
            (
                "Fnet",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ
            ),
            mesh,
            dimensionedScalar(dimForce, 0),
            zeroGradientFvPatchScalarField::typeName
        );
        volFnet.primitiveFieldRef() = Fnet;
        volFnet.correctBoundaryConditions();
        volFnet.write();
    }

    ejectionModel::correct();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace surfaceFilmModels
} // End namespace regionModels
} // End namespace Foam

// ************************************************************************* //
