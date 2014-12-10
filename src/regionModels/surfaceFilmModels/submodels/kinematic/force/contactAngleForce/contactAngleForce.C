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

#include "contactAngleForce.H"
#include "addToRunTimeSelectionTable.H"
#include "fvcGrad.H"
#include "unitConversion.H"
#include "fvPatchField.H"
#include "patchDist.H"

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
            label patchI = iter.key();
            Info<< "            " << pbm[patchI].name() << endl;
        }

        patchDist dist(owner_.regionMesh(), patchIDs);

        mask_ = pos(dist - dimensionedScalar("dLim", dimLength, dLim));
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
            dimensionedVector("zero", dimForce/dimArea, vector::zero)
        )
    );

    vectorField& force = tForce().internalField();

    const labelUList& own = owner_.regionMesh().owner();
    const labelUList& nbr = owner_.regionMesh().neighbour();

    const scalarField& magSf = owner_.magSf();

    const volScalarField& alpha = owner_.alpha();
    const volScalarField& sigma = owner_.sigma();

    volVectorField gradAlpha(fvc::grad(alpha));

    forAll(nbr, faceI)
    {
        const label cellO = own[faceI];
        const label cellN = nbr[faceI];

        label cellI = -1;
        if ((alpha[cellO] > 0.5) && (alpha[cellN] < 0.5))
        {
            cellI = cellO;
        }
        else if ((alpha[cellO] < 0.5) && (alpha[cellN] > 0.5))
        {
            cellI = cellN;
        }

        if (cellI != -1 && mask_[cellI] > 0.5)
        {
            const scalar invDx = owner_.regionMesh().deltaCoeffs()[faceI];
            const vector n =
                gradAlpha[cellI]/(mag(gradAlpha[cellI]) + ROOTVSMALL);
            scalar theta = cos(degToRad(distribution_->sample()));
            force[cellI] += Ccf_*n*sigma[cellI]*(1.0 - theta)/invDx;
        }
    }

    forAll(alpha.boundaryField(), patchI)
    {
        if (!owner().isCoupledPatch(patchI))
        {
            const fvPatchField<scalar>& alphaf = alpha.boundaryField()[patchI];
            const fvPatchField<scalar>& maskf = mask_.boundaryField()[patchI];
            const scalarField& invDx = alphaf.patch().deltaCoeffs();
            const labelUList& faceCells = alphaf.patch().faceCells();

            forAll(alphaf, faceI)
            {
                if (maskf[faceI] > 0.5)
                {
                    label cellO = faceCells[faceI];

                    if ((alpha[cellO] > 0.5) && (alphaf[faceI] < 0.5))
                    {
                        const vector n =
                            gradAlpha[cellO]
                           /(mag(gradAlpha[cellO]) + ROOTVSMALL);
                        scalar theta = cos(degToRad(distribution_->sample()));
                        force[cellO] +=
                            Ccf_*n*sigma[cellO]*(1.0 - theta)/invDx[faceI];
                    }
                }
            }
        }
    }

    force /= magSf;

    if (owner_.regionMesh().time().outputTime())
    {
        tForce().write();
    }

    tmp<fvVectorMatrix>
        tfvm(new fvVectorMatrix(U, dimForce/dimArea*dimVolume));

    tfvm() += tForce;

    return tfvm;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace surfaceFilmModels
} // End namespace regionModels
} // End namespace Foam

// ************************************************************************* //
