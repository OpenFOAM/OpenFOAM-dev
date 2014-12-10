/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

#include "Wallis.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace compressibilityModels
    {
        defineTypeNameAndDebug(Wallis, 0);
        addToRunTimeSelectionTable
        (
            barotropicCompressibilityModel,
            Wallis,
            dictionary
        );
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::compressibilityModels::Wallis::Wallis
(
    const dictionary& compressibilityProperties,
    const volScalarField& gamma,
    const word& psiName
)
:
    barotropicCompressibilityModel(compressibilityProperties, gamma, psiName),
    psiv_(compressibilityProperties_.lookup("psiv")),
    psil_(compressibilityProperties_.lookup("psil")),
    rhovSat_(compressibilityProperties_.lookup("rhovSat")),
    rholSat_(compressibilityProperties_.lookup("rholSat"))
{
    correct();
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::compressibilityModels::Wallis::correct()
{
    psi_ =
        (gamma_*rhovSat_ + (scalar(1) - gamma_)*rholSat_)
       *(gamma_*psiv_/rhovSat_ + (scalar(1) - gamma_)*psil_/rholSat_);
}


bool Foam::compressibilityModels::Wallis::read
(
    const dictionary& compressibilityProperties
)
{
    barotropicCompressibilityModel::read(compressibilityProperties);

    compressibilityProperties_.lookup("psiv") >> psiv_;
    compressibilityProperties_.lookup("psil") >> psil_;
    compressibilityProperties_.lookup("rhovSat") >> rhovSat_;
    compressibilityProperties_.lookup("rholSat") >> rholSat_;

    return true;
}


// ************************************************************************* //
