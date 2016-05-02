/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015-2016 OpenFOAM Foundation
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

#include "yPlus.H"
#include "nutWallFunctionFvPatchScalarField.H"
#include "nearWallDist.H"
#include "wallFvPatch.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class TurbulenceModel>
void Foam::functionObjects::yPlus::calcYPlus
(
    const TurbulenceModel& turbulenceModel,
    const fvMesh& mesh,
    volScalarField& yPlus
)
{
    volScalarField::Boundary d = nearWallDist(mesh).y();

    const volScalarField::Boundary nutBf =
        turbulenceModel.nut()().boundaryField();

    const volScalarField::Boundary nuEffBf =
        turbulenceModel.nuEff()().boundaryField();

    const volScalarField::Boundary nuBf =
        turbulenceModel.nu()().boundaryField();

    const fvPatchList& patches = mesh.boundary();

    volScalarField::Boundary& yPlusBf =
        yPlus.boundaryFieldRef();

    forAll(patches, patchi)
    {
        const fvPatch& patch = patches[patchi];

        if (isA<nutWallFunctionFvPatchScalarField>(nutBf[patchi]))
        {
            const nutWallFunctionFvPatchScalarField& nutPf =
                dynamic_cast<const nutWallFunctionFvPatchScalarField&>
                (
                    nutBf[patchi]
                );

            yPlusBf[patchi] = nutPf.yPlus();
            const scalarField& yPlusp = yPlusBf[patchi];

            const scalar minYplus = gMin(yPlusp);
            const scalar maxYplus = gMax(yPlusp);
            const scalar avgYplus = gAverage(yPlusp);

            if (Pstream::master())
            {
                if (log_) Info
                    << "    patch " << patch.name()
                    << " y+ : min = " << minYplus << ", max = " << maxYplus
                    << ", average = " << avgYplus << nl;

                writeTime(file());
                file()
                    << token::TAB << patch.name()
                    << token::TAB << minYplus
                    << token::TAB << maxYplus
                    << token::TAB << avgYplus
                    << endl;
            }
        }
        else if (isA<wallFvPatch>(patch))
        {
            yPlusBf[patchi] =
                d[patchi]
               *sqrt
                (
                    nuEffBf[patchi]
                   *mag(turbulenceModel.U().boundaryField()[patchi].snGrad())
                )/nuBf[patchi];
            const scalarField& yPlusp = yPlusBf[patchi];

            const scalar minYplus = gMin(yPlusp);
            const scalar maxYplus = gMax(yPlusp);
            const scalar avgYplus = gAverage(yPlusp);

            if (Pstream::master())
            {
                if (log_) Info
                    << "    patch " << patch.name()
                    << " y+ : min = " << minYplus << ", max = " << maxYplus
                    << ", average = " << avgYplus << nl;

                writeTime(file());
                file()
                    << token::TAB << patch.name()
                    << token::TAB << minYplus
                    << token::TAB << maxYplus
                    << token::TAB << avgYplus
                    << endl;
            }
        }
    }
}


// ************************************************************************* //
