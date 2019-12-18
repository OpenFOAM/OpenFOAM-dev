/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
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

#include "Chung.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace compressibilityModels
    {
        defineTypeNameAndDebug(Chung, 0);
        addToRunTimeSelectionTable
        (
            barotropicCompressibilityModel,
            Chung,
            dictionary
        );
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::compressibilityModels::Chung::Chung
(
    const dictionary& compressibilityProperties,
    const volScalarField& gamma,
    const word& psiName
)
:
    barotropicCompressibilityModel(compressibilityProperties, gamma, psiName),
    pSat_
    (
        "pSat",
        dimPressure,
        compressibilityProperties_.lookup("pSat")
    ),
    psiv_
    (
        "psiv",
        dimCompressibility,
        compressibilityProperties_.lookup("psiv")
    ),
    psil_
    (
        "psil",
        dimCompressibility,
        compressibilityProperties_.lookup("psil")
    ),
    rholSat_
    (
        "rholSat",
        dimDensity,
        compressibilityProperties_.lookup("rholSat")
    )
{
    correct();
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::compressibilityModels::Chung::correct()
{
    volScalarField sfa
    (
        sqrt
        (
            pSat_
           /((scalar(1) - gamma_)*pSat_ + gamma_*rholSat_/psil_)
        )
    );

    psi_ = sqr
    (
        ((scalar(1) - gamma_)/sqrt(psiv_) + gamma_*sfa/sqrt(psil_))
       *sqrt(psiv_*psil_)/sfa
    );
}


bool Foam::compressibilityModels::Chung::read
(
    const dictionary& compressibilityProperties
)
{
    barotropicCompressibilityModel::read(compressibilityProperties);

    compressibilityProperties_.lookup("pSat") >> pSat_;
    compressibilityProperties_.lookup("psiv") >> psiv_;
    compressibilityProperties_.lookup("psil") >> psil_;
    compressibilityProperties_.lookup("rholSat") >> rholSat_;

    return true;
}


// ************************************************************************* //
