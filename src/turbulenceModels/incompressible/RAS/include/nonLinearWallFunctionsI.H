/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
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

Global
    nonLinearwallFunctions

Description
    Calculate wall generation and dissipation from wall-functions
    for non-linear models.

\*---------------------------------------------------------------------------*/

{
    labelList cellBoundaryFaceCount(epsilon_.size(), 0);

    scalar yPlusLam = nutkWallFunctionFvPatchScalarField::yPlusLam
    (
        kappa_.value(),
        E_.value()
    );

    const fvPatchList& patches = mesh_.boundary();

    //- Initialise the near-wall G and epsilon fields to zero
    forAll(patches, patchi)
    {
        const fvPatch& curPatch = patches[patchi];

        if (isA<wallFvPatch>(curPatch))
        {
            forAll(curPatch, facei)
            {
                label faceCelli = curPatch.faceCells()[facei];

                epsilon_[faceCelli] = 0.0;
                G[faceCelli] = 0.0;
            }
        }
    }

    const volScalarField nuLam(this->nu());

    //- Accumulate the wall face contributions to epsilon and G
    //  Increment cellBoundaryFaceCount for each face for averaging
    forAll(patches, patchi)
    {
        const fvPatch& curPatch = patches[patchi];

        if (isA<wallFvPatch>(curPatch))
        {
            #include "checkPatchFieldTypes.H"

            const scalarField& nuw = nuLam.boundaryField()[patchi];
            const scalarField& nutw = nut_.boundaryField()[patchi];

            const scalarField magFaceGradU
            (
                mag(U_.boundaryField()[patchi].snGrad())
            );

            forAll(curPatch, facei)
            {
                label faceCelli = curPatch.faceCells()[facei];

                //- using local Cmu !
                scalar Cmu25 = pow025(Cmu_[faceCelli]);
                scalar Cmu75 = pow(Cmu_[faceCelli], 0.75);

                scalar yPlus =
                    Cmu25*y_[patchi][facei]
                    *sqrt(k_[faceCelli])
                    /nuw[facei];

                // For corner cells (with two boundary or more faces),
                // epsilon and G in the near-wall cell are calculated
                // as an average

                cellBoundaryFaceCount[faceCelli]++;

                epsilon_[faceCelli] +=
                     Cmu75*pow(k_[faceCelli], 1.5)
                    /(kappa_.value()*y_[patchi][facei]);

                if (yPlus > yPlusLam)
                {
                    G[faceCelli] +=
                        (nutw[facei] + nuw[facei])
                        *magFaceGradU[facei]
                        *Cmu25*sqrt(k_[faceCelli])
                        /(kappa_.value()*y_[patchi][facei])
                      - (nonlinearStress_[faceCelli] && gradU[faceCelli]);
                }
            }
        }
    }

    // Perform the averaging

    forAll(patches, patchi)
    {
        const fvPatch& curPatch = patches[patchi];

        if (isA<wallFvPatch>(curPatch))
        {
            forAll(curPatch, facei)
            {
                label faceCelli = curPatch.faceCells()[facei];

                epsilon_[faceCelli] /= cellBoundaryFaceCount[faceCelli];
                G[faceCelli] /= cellBoundaryFaceCount[faceCelli];
            }
        }
    }
}


// ************************************************************************* //
