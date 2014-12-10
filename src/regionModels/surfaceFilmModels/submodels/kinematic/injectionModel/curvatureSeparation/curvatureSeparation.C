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

#include "curvatureSeparation.H"
#include "addToRunTimeSelectionTable.H"
#include "fvMesh.H"
#include "Time.H"
#include "volFields.H"
#include "kinematicSingleLayer.H"
#include "surfaceInterpolate.H"
#include "fvcDiv.H"
#include "fvcGrad.H"
#include "stringListOps.H"
#include "cyclicPolyPatch.H"

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
    injectionModel,
    curvatureSeparation,
    dictionary
);

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

tmp<volScalarField> curvatureSeparation::calcInvR1
(
    const volVectorField& U
) const
{
    // method 1
/*
    tmp<volScalarField> tinvR1
    (
        new volScalarField("invR1", fvc::div(owner().nHat()))
    );
*/

    // method 2
    dimensionedScalar smallU("smallU", dimVelocity, ROOTVSMALL);
    volVectorField UHat(U/(mag(U) + smallU));
    tmp<volScalarField> tinvR1
    (
        new volScalarField("invR1", UHat & (UHat & gradNHat_))
    );


    scalarField& invR1 = tinvR1().internalField();

    // apply defined patch radii
    const scalar rMin = 1e-6;
    const fvMesh& mesh = owner().regionMesh();
    const polyBoundaryMesh& pbm = mesh.boundaryMesh();
    forAll(definedPatchRadii_, i)
    {
        label patchI = definedPatchRadii_[i].first();
        scalar definedInvR1 = 1.0/max(rMin, definedPatchRadii_[i].second());
        UIndirectList<scalar>(invR1, pbm[patchI].faceCells()) = definedInvR1;
    }

    // filter out large radii
    const scalar rMax = 1e6;
    forAll(invR1, i)
    {
        if (mag(invR1[i]) < 1/rMax)
        {
            invR1[i] = -1.0;
        }
    }

    if (debug && mesh.time().outputTime())
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
    const fvMesh& mesh = owner().regionMesh();
    const vectorField nf(mesh.Sf()/mesh.magSf());
    const unallocLabelList& own = mesh.owner();
    const unallocLabelList& nbr = mesh.neighbour();

    scalarField phiMax(mesh.nCells(), -GREAT);
    scalarField cosAngle(mesh.nCells(), 0.0);
    forAll(nbr, faceI)
    {
        label cellO = own[faceI];
        label cellN = nbr[faceI];

        if (phi[faceI] > phiMax[cellO])
        {
            phiMax[cellO] = phi[faceI];
            cosAngle[cellO] = -gHat_ & nf[faceI];
        }
        if (-phi[faceI] > phiMax[cellN])
        {
            phiMax[cellN] = -phi[faceI];
            cosAngle[cellN] = -gHat_ & -nf[faceI];
        }
    }

    forAll(phi.boundaryField(), patchI)
    {
        const fvsPatchScalarField& phip = phi.boundaryField()[patchI];
        const fvPatch& pp = phip.patch();
        const labelList& faceCells = pp.faceCells();
        const vectorField nf(pp.nf());
        forAll(phip, i)
        {
            label cellI = faceCells[i];
            if (phip[i] > phiMax[cellI])
            {
                phiMax[cellI] = phip[i];
                cosAngle[cellI] = -gHat_ & nf[i];
            }
        }
    }
/*
    // correction for cyclics - use cyclic pairs' face normal instead of
    // local face normal
    const fvBoundaryMesh& pbm = mesh.boundary();
    forAll(phi.boundaryField(), patchI)
    {
        if (isA<cyclicPolyPatch>(pbm[patchI]))
        {
            const scalarField& phip = phi.boundaryField()[patchI];
            const vectorField nf(pbm[patchI].nf());
            const labelList& faceCells = pbm[patchI].faceCells();
            const label sizeBy2 = pbm[patchI].size()/2;

            for (label face0=0; face0<sizeBy2; face0++)
            {
                label face1 = face0 + sizeBy2;
                label cell0 = faceCells[face0];
                label cell1 = faceCells[face1];

                // flux leaving half 0, entering half 1
                if (phip[face0] > phiMax[cell0])
                {
                    phiMax[cell0] = phip[face0];
                    cosAngle[cell0] = -gHat_ & -nf[face1];
                }

                // flux leaving half 1, entering half 0
                if (-phip[face1] > phiMax[cell1])
                {
                    phiMax[cell1] = -phip[face1];
                    cosAngle[cell1] = -gHat_ & nf[face0];
                }
            }
        }
    }
*/
    // checks
    if (debug && mesh.time().outputTime())
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
            dimensionedScalar("zero", dimless, 0.0),
            zeroGradientFvPatchScalarField::typeName
        );
        volCosAngle.internalField() = cosAngle;
        volCosAngle.correctBoundaryConditions();
        volCosAngle.write();
    }

    return max(min(cosAngle, scalar(1.0)), scalar(-1.0));
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

curvatureSeparation::curvatureSeparation
(
    surfaceFilmModel& owner,
    const dictionary& dict
)
:
    injectionModel(type(), owner, dict),
    gradNHat_(fvc::grad(owner.nHat())),
    deltaByR1Min_(coeffDict_.lookupOrDefault<scalar>("deltaByR1Min", 0.0)),
    definedPatchRadii_(),
    magG_(mag(owner.g().value())),
    gHat_(vector::zero)
{
    if (magG_ < ROOTVSMALL)
    {
        FatalErrorIn
        (
            "curvatureSeparation::curvatureSeparation"
            "("
                "const surfaceFilmModel&, "
                "const dictionary&"
            ")"
        )
            << "Acceleration due to gravity must be non-zero"
            << exit(FatalError);
    }

    gHat_ = owner.g().value()/magG_;

    List<Tuple2<word, scalar> > prIn(coeffDict_.lookup("definedPatchRadii"));
    const wordList& allPatchNames = owner.regionMesh().boundaryMesh().names();

    DynamicList<Tuple2<label, scalar> > prData(allPatchNames.size());

    labelHashSet uniquePatchIDs;

    forAllReverse(prIn, i)
    {
        labelList patchIDs = findStrings(prIn[i].first(), allPatchNames);
        forAll(patchIDs, j)
        {
            const label patchI = patchIDs[j];

            if (!uniquePatchIDs.found(patchI))
            {
                const scalar radius = prIn[i].second();
                prData.append(Tuple2<label, scalar>(patchI, radius));

                uniquePatchIDs.insert(patchI);
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
    scalarField& massToInject,
    scalarField& diameterToInject
)
{
    const kinematicSingleLayer& film =
        refCast<const kinematicSingleLayer>(this->owner());
    const fvMesh& mesh = film.regionMesh();

    const volScalarField& delta = film.delta();
    const volVectorField& U = film.U();
    const surfaceScalarField& phi = film.phi();
    const volScalarField& rho = film.rho();
    const scalarField magSqrU(magSqr(film.U()));
    const volScalarField& sigma = film.sigma();

    const scalarField invR1(calcInvR1(U));
    const scalarField cosAngle(calcCosAngle(phi));

    // calculate force balance
    const scalar Fthreshold = 1e-10;
    scalarField Fnet(mesh.nCells(), 0.0);
    scalarField separated(mesh.nCells(), 0.0);
    forAll(invR1, i)
    {
        if ((invR1[i] > 0) && (delta[i]*invR1[i] > deltaByR1Min_))
        {
            scalar R1 = 1.0/(invR1[i] + ROOTVSMALL);
            scalar R2 = R1 + delta[i];

            // inertial force
            scalar Fi = -delta[i]*rho[i]*magSqrU[i]*72.0/60.0*invR1[i];

            // body force
            scalar Fb =
              - 0.5*rho[i]*magG_*invR1[i]*(sqr(R1) - sqr(R2))*cosAngle[i];

            // surface force
            scalar Fs = sigma[i]/R2;

            Fnet[i] = Fi + Fb + Fs;

            if (Fnet[i] + Fthreshold < 0)
            {
                separated[i] = 1.0;
            }
        }
    }

    // inject all available mass
    massToInject = separated*availableMass;
    diameterToInject = separated*delta;
    availableMass -= separated*availableMass;

    addToInjectedMass(sum(separated*availableMass));

    if (debug && mesh.time().outputTime())
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
            dimensionedScalar("zero", dimForce, 0.0),
            zeroGradientFvPatchScalarField::typeName
        );
        volFnet.internalField() = Fnet;
        volFnet.correctBoundaryConditions();
        volFnet.write();
    }

    injectionModel::correct();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace surfaceFilmModels
} // End namespace regionModels
} // End namespace Foam

// ************************************************************************* //
