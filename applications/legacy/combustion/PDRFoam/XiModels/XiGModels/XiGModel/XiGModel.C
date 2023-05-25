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

#include "XiGModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(XiGModel, 0);
    defineRunTimeSelectionTable(XiGModel, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::XiGModel::XiGModel
(
    const dictionary& XiGProperties,
    const psiuMulticomponentThermo& thermo,
    const compressible::RASModel& turbulence,
    const volScalarField& Su
)
:
    XiGModelCoeffs_
    (
        XiGProperties.subDict
        (
            word(XiGProperties.lookup("XiGModel")) + "Coeffs"
        )
    ),
    thermo_(thermo),
    turbulence_(turbulence),
    Su_(Su)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::XiGModel::~XiGModel()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::XiGModel::read(const dictionary& XiGProperties)
{
    XiGModelCoeffs_ = XiGProperties.optionalSubDict(type() + "Coeffs");

    return true;
}


// ************************************************************************* //
