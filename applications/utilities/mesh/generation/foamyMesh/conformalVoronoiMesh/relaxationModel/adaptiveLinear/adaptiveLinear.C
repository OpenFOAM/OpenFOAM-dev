/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2012-2019 OpenFOAM Foundation
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

#include "adaptiveLinear.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(adaptiveLinear, 0);
addToRunTimeSelectionTable(relaxationModel, adaptiveLinear, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

adaptiveLinear::adaptiveLinear
(
    const dictionary& relaxationDict,
    const Time& runTime
)
:
    relaxationModel(typeName, relaxationDict, runTime),
    relaxationStart_(coeffDict().lookup<scalar>("relaxationStart")),
    relaxationEnd_(coeffDict().lookup<scalar>("relaxationEnd")),
    lastTimeValue_(runTime_.time().timeOutputValue()),
    relaxation_(relaxationStart_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

scalar adaptiveLinear::relaxation()
{
    if (runTime_.time().timeOutputValue() > lastTimeValue_)
    {
        scalar currentRelaxation = relaxation_;

        relaxation_ -=
            (relaxation_ - relaxationEnd_)
           /(
                (
                    runTime_.time().endTime().value()
                  - runTime_.time().timeOutputValue()
                )
               /(runTime_.time().timeOutputValue() - lastTimeValue_)
              + 1
            );

        lastTimeValue_ = runTime_.time().timeOutputValue();

        return currentRelaxation;
    }

    return relaxation_;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
