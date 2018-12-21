/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

void contactAngleForce::initialise()
{
    const wordReList zeroForcePatches
    (
        coeffDict_.lookupOrDefault<wordReList>("zeroForcePatches", wordReList())
    );

    if (zeroForcePatches.size())
    {
        const polyBoundaryMesh& pbm = filmModel_.regionMesh().boundaryMesh();
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
        patchDistMethods::meshWave dist(filmModel_.regionMesh(), patchIDs);
        volScalarField y
        (
            IOobject
            (
                "y",
                filmModel_.regionMesh().time().timeName(),
                filmModel_.regionMesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            filmModel_.regionMesh(),
            dimensionedScalar(dimLength, great)
        );
        dist.correct(y);

        mask_ = pos0(y - dimensionedScalar(dimLength, dLim));
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

contactAngleForce::contactAngleForce
(
    const word& typeName,
    surfaceFilmRegionModel& film,
    const dictionary& dict
)
:
    force(typeName, film, dict),
    Ccf_(readScalar(coeffDict_.lookup("Ccf"))),
    mask_
    (
        IOobject
        (
            typeName + ":contactForceMask",
            filmModel_.time().timeName(),
            filmModel_.regionMesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        filmModel_.regionMesh(),
        dimensionedScalar(dimless, 1.0)
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
        volVectorField::New
        (
            typeName + ":contactForce",
            filmModel_.regionMesh(),
            dimensionedVector(dimForce/dimArea, Zero)
        )
    );

    vectorField& force = tForce.ref().primitiveFieldRef();

    const labelUList& own = filmModel_.regionMesh().owner();
    const labelUList& nbr = filmModel_.regionMesh().neighbour();

    const scalarField& magSf = filmModel_.magSf();

    const volScalarField& alpha = filmModel_.alpha();
    const volScalarField& sigma = filmModel_.sigma();

    const tmp<volScalarField> ttheta = theta();
    const volScalarField& theta = ttheta();

    const volVectorField gradAlpha(fvc::grad(alpha));

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
            const scalar invDx = filmModel_.regionMesh().deltaCoeffs()[facei];
            const vector n =
                gradAlpha[celli]/(mag(gradAlpha[celli]) + rootVSmall);
            const scalar cosTheta = cos(degToRad(theta[celli]));
            force[celli] += Ccf_*n*sigma[celli]*(1 - cosTheta)/invDx;
        }
    }

    forAll(alpha.boundaryField(), patchi)
    {
        if (!filmModel_.isCoupledPatch(patchi))
        {
            const fvPatchField<scalar>& alphaPf = alpha.boundaryField()[patchi];
            const fvPatchField<scalar>& maskPf = mask_.boundaryField()[patchi];
            const fvPatchField<scalar>& sigmaPf = sigma.boundaryField()[patchi];
            const fvPatchField<scalar>& thetaPf = theta.boundaryField()[patchi];
            const scalarField& invDx = alphaPf.patch().deltaCoeffs();
            const labelUList& faceCells = alphaPf.patch().faceCells();

            forAll(alphaPf, facei)
            {
                if (maskPf[facei] > 0.5)
                {
                    label cellO = faceCells[facei];

                    if ((alpha[cellO] > 0.5) && (alphaPf[facei] < 0.5))
                    {
                        const vector n =
                            gradAlpha[cellO]
                           /(mag(gradAlpha[cellO]) + rootVSmall);
                        const scalar cosTheta = cos(degToRad(thetaPf[facei]));
                        force[cellO] +=
                            Ccf_*n*sigmaPf[facei]*(1 - cosTheta)/invDx[facei];
                    }
                }
            }
        }
    }

    force /= magSf;

    if (filmModel_.regionMesh().time().writeTime())
    {
        tForce().write();
    }

    tmp<fvVectorMatrix> tfvm
    (
        new fvVectorMatrix(U, dimForce/dimArea*dimVolume)
    );

    tfvm.ref() += tForce;

    return tfvm;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace surfaceFilmModels
} // End namespace regionModels
} // End namespace Foam

// ************************************************************************* //
