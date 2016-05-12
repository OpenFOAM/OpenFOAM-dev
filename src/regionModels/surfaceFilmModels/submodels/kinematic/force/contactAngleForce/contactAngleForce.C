/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
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

#include "contactAngleForce.H"
#include "addToRunTimeSelectionTable.H"
#include "fvcGrad.H"
#include "unitConversion.H"
#include "fvPatchField.H"
#include "meshWavePatchDistMethod.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace surfaceFilmModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(contactAngleForce, 0);
addToRunTimeSelectionTable(force, contactAngleForce, dictionary);

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

void contactAngleForce::initialise()
{
    const wordReList zeroForcePatches(coeffDict_.lookup("zeroForcePatches"));

    if (zeroForcePatches.size())
    {
        const polyBoundaryMesh& pbm = owner_.regionMesh().boundaryMesh();
        scalar dLim = readScalar(coeffDict_.lookup("zeroForceDistance"));

        Info<< "        Assigning zero contact force within " << dLim
            << " of patches:" << endl;

        labelHashSet patchIDs = pbm.patchSet(zeroForcePatches);

        forAllConstIter(labelHashSet, patchIDs, iter)
        {
            label patchi = iter.key();
            Info<< "            " << pbm[patchi].name() << endl;
        }

        // Temporary implementation until run-time selection covers this case
        patchDistMethods::meshWave dist(owner_.regionMesh(), patchIDs);
        volScalarField y
        (
            IOobject
            (
                "y",
                owner_.regionMesh().time().timeName(),
                owner_.regionMesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            owner_.regionMesh(),
            dimensionedScalar("y", dimLength, GREAT)
        );
        dist.correct(y);

        mask_ = pos(y - dimensionedScalar("dLim", dimLength, dLim));
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

contactAngleForce::contactAngleForce
(
    surfaceFilmModel& owner,
    const dictionary& dict
)
:
    force(typeName, owner, dict),
    Ccf_(readScalar(coeffDict_.lookup("Ccf"))),
    rndGen_(label(0), -1),
    distribution_
    (
        distributionModels::distributionModel::New
        (
            coeffDict_.subDict("contactAngleDistribution"),
            rndGen_
        )
    ),
    mask_
    (
        IOobject
        (
            typeName + ":contactForceMask",
            owner_.time().timeName(),
            owner_.regionMesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        owner_.regionMesh(),
        dimensionedScalar("mask", dimless, 1.0)
    )
{
    initialise();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

contactAngleForce::~contactAngleForce()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

tmp<fvVectorMatrix> contactAngleForce::correct(volVectorField& U)
{
    tmp<volVectorField> tForce
    (
        new volVectorField
        (
            IOobject
            (
                typeName + ":contactForce",
                owner_.time().timeName(),
                owner_.regionMesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            owner_.regionMesh(),
            dimensionedVector("zero", dimForce/dimArea, Zero)
        )
    );

    vectorField& force = tForce.ref().primitiveFieldRef();

    const labelUList& own = owner_.regionMesh().owner();
    const labelUList& nbr = owner_.regionMesh().neighbour();

    const scalarField& magSf = owner_.magSf();

    const volScalarField& alpha = owner_.alpha();
    const volScalarField& sigma = owner_.sigma();

    volVectorField gradAlpha(fvc::grad(alpha));

    forAll(nbr, facei)
    {
        const label cellO = own[facei];
        const label cellN = nbr[facei];

        label celli = -1;
        if ((alpha[cellO] > 0.5) && (alpha[cellN] < 0.5))
        {
            celli = cellO;
        }
        else if ((alpha[cellO] < 0.5) && (alpha[cellN] > 0.5))
        {
            celli = cellN;
        }

        if (celli != -1 && mask_[celli] > 0.5)
        {
            const scalar invDx = owner_.regionMesh().deltaCoeffs()[facei];
            const vector n =
                gradAlpha[celli]/(mag(gradAlpha[celli]) + ROOTVSMALL);
            scalar theta = cos(degToRad(distribution_->sample()));
            force[celli] += Ccf_*n*sigma[celli]*(1.0 - theta)/invDx;
        }
    }

    forAll(alpha.boundaryField(), patchi)
    {
        if (!owner().isCoupledPatch(patchi))
        {
            const fvPatchField<scalar>& alphaf = alpha.boundaryField()[patchi];
            const fvPatchField<scalar>& maskf = mask_.boundaryField()[patchi];
            const scalarField& invDx = alphaf.patch().deltaCoeffs();
            const labelUList& faceCells = alphaf.patch().faceCells();

            forAll(alphaf, facei)
            {
                if (maskf[facei] > 0.5)
                {
                    label cellO = faceCells[facei];

                    if ((alpha[cellO] > 0.5) && (alphaf[facei] < 0.5))
                    {
                        const vector n =
                            gradAlpha[cellO]
                           /(mag(gradAlpha[cellO]) + ROOTVSMALL);
                        scalar theta = cos(degToRad(distribution_->sample()));
                        force[cellO] +=
                            Ccf_*n*sigma[cellO]*(1.0 - theta)/invDx[facei];
                    }
                }
            }
        }
    }

    force /= magSf;

    if (owner_.regionMesh().time().writeTime())
    {
        tForce().write();
    }

    tmp<fvVectorMatrix>
        tfvm(new fvVectorMatrix(U, dimForce/dimArea*dimVolume));

    tfvm.ref() += tForce;

    return tfvm;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace surfaceFilmModels
} // End namespace regionModels
} // End namespace Foam

// ************************************************************************* //
