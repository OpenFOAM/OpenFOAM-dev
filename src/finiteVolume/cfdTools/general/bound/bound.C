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

#include "bound.H"
#include "fvcAverage.H"

// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

bool Foam::bound(volScalarField& vsf, const dimensionedScalar& min)
{
    const scalar minVsf = Foam::min(vsf).value();

    if (minVsf < min.value())
    {
        scalarField& isf = vsf.primitiveFieldRef();

        Info<< "bounding " << vsf.name()
            << ", min: " << minVsf
            << " max: " << max(vsf).value()
            << " average: " << gAverage(isf)
            << endl;

        isf = max
        (
            max
            (
                isf,
                fvc::average(max(vsf, min))().primitiveField()
               *pos0(-isf)
            ),
            min.value()
        );

        volScalarField::Boundary& bsf = vsf.boundaryFieldRef();
        forAll(bsf, patchi)
        {
            bsf[patchi] == max
            (
                max
                (
                    bsf[patchi],
                    bsf[patchi].patchInternalField()
                   *pos0(-bsf[patchi])
                ),
                min.value()
            );
        }

        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
