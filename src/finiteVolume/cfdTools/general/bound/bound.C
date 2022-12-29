/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2022 OpenFOAM Foundation
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
#include "volFields.H"
#include "fvcAverage.H"

// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

Foam::volScalarField&
Foam::bound(volScalarField& vsf, const dimensionedScalar& lowerBound)
{
    const scalar minVsf = min(vsf).value();

    if (minVsf < lowerBound.value())
    {
        Info<< "bounding " << vsf.name()
            << ", min: " << minVsf
            << " max: " << max(vsf).value()
            << " average: " << gAverage(vsf.primitiveField())
            << endl;

        vsf.primitiveFieldRef() = max
        (
            max
            (
                vsf.primitiveField(),
                fvc::average(max(vsf, lowerBound))().primitiveField()
              * pos0(-vsf.primitiveField())
            ),
            lowerBound.value()
        );

        vsf.boundaryFieldRef() = max(vsf.boundaryField(), lowerBound.value());
    }

    return vsf;
}


// ************************************************************************* //
