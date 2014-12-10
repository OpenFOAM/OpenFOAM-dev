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
    wallNonlinearViscosity

Description
    Calculate wall viscosity for non-linear models

\*---------------------------------------------------------------------------*/

{
    const fvPatchList& patches = mesh_.boundary();

    const scalar yPlusLam = nutkWallFunctionFvPatchScalarField::yPlusLam
    (
        kappa_.value(),
        E_.value()
    );

    const volScalarField nuLam(this->nu());

    forAll(patches, patchi)
    {
        const fvPatch& curPatch = patches[patchi];

        if (isA<wallFvPatch>(curPatch))
        {
            const scalarField& nuw = nuLam.boundaryField()[patchi];
            scalarField& nutw = nut_.boundaryField()[patchi];

            forAll(curPatch, facei)
            {
                label faceCelli = curPatch.faceCells()[facei];

                //- Using local Cmu
                scalar Cmu25 = pow025(Cmu_[faceCelli]);

                scalar yPlus =
                    Cmu25*y_[patchi][facei]*sqrt(k_[faceCelli])/nuw[facei];

                if (yPlus > yPlusLam)
                {
                    nutw[facei] =
                        nuw[facei]
                       *(yPlus*kappa_.value()/log(E_.value()*yPlus) - 1);
                }
                else
                {
                    nutw[facei] = 0.0;
                }
            }
        }
    }
}


// ************************************************************************* //
