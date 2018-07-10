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

InClass
    barotropicCompressibilityModel

\*---------------------------------------------------------------------------*/

#include "barotropicCompressibilityModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(barotropicCompressibilityModel, 0);
    defineRunTimeSelectionTable(barotropicCompressibilityModel, dictionary);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::barotropicCompressibilityModel::barotropicCompressibilityModel
(
    const dictionary& compressibilityProperties,
    const volScalarField& gamma,
    const word& psiName
)
:
    compressibilityProperties_(compressibilityProperties),
    psi_
    (
        IOobject
        (
            psiName,
            gamma.mesh().time().timeName(),
            gamma.mesh()
        ),
        gamma.mesh(),
        dimensionedScalar(psiName, dimensionSet(0, -2, 2, 0, 0), 0)
    ),
    gamma_(gamma)
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::barotropicCompressibilityModel::read
(
    const dictionary& compressibilityProperties
)
{
    compressibilityProperties_ = compressibilityProperties;

    return true;
}


// ************************************************************************* //
