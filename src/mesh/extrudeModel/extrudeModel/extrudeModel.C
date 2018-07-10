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

\*---------------------------------------------------------------------------*/

#include "extrudeModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(extrudeModel, 0);
    defineRunTimeSelectionTable(extrudeModel, dictionary);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::extrudeModel::extrudeModel
(
    const word& modelType,
    const dictionary& dict
)
:
    nLayers_(dict.lookupOrDefault<label>("nLayers", 1)),
    expansionRatio_(dict.lookupOrDefault<scalar>("expansionRatio", 1)),
    dict_(dict),
    coeffDict_(dict.optionalSubDict(modelType + "Coeffs"))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::extrudeModel::~extrudeModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::label Foam::extrudeModel::nLayers() const
{
    return nLayers_;
}


Foam::scalar Foam::extrudeModel::expansionRatio() const
{
    return expansionRatio_;
}


Foam::scalar Foam::extrudeModel::sumThickness(const label layer) const
{
    // 1+r+r^2+ .. +r^(n-1) = (1-r^n)/(1-r)

    if (mag(1.0-expansionRatio_) < small)
    {
        return scalar(layer)/nLayers_;
    }
    else
    {
        return
            (1.0-pow(expansionRatio_, layer))
          / (1.0-pow(expansionRatio_, nLayers_));
    }
}


// ************************************************************************* //
