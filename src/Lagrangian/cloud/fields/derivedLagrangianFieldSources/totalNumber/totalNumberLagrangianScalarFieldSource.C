/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2025 OpenFOAM Foundation
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

#include "totalNumberLagrangianScalarFieldSource.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::totalNumberLagrangianScalarFieldSource::
totalNumberLagrangianScalarFieldSource
(
    const regIOobject& iIo,
    const dictionary& dict
)
:
    uniformSizeNumberLagrangianScalarFieldSource(iIo, dict),
    haveVolume_(dict.found("volume")),
    volumeOrMass_
    (
        haveVolume_
      ? dict.lookup<scalar>("volume", dimVolume)
      : dict.found("mass")
      ? dict.lookup<scalar>("mass", dimMass)
      : NaN
    )
{
    bool haveTotalMass = dict.found("mass");

    if (haveVolume_ == haveTotalMass)
    {
        FatalIOErrorInFunction(dict)
            << "keywords volume and mass both "
            << (haveVolume_ ? "" : "un") << "defined in "
            << "dictionary " << dict.name()
            << exit(FatalIOError);
    }
}


Foam::totalNumberLagrangianScalarFieldSource::
totalNumberLagrangianScalarFieldSource
(
    const totalNumberLagrangianScalarFieldSource& field,
    const regIOobject& iIo
)
:
    uniformSizeNumberLagrangianScalarFieldSource(field, iIo),
    haveVolume_(field.haveVolume_),
    volumeOrMass_(field.volumeOrMass_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::totalNumberLagrangianScalarFieldSource::
~totalNumberLagrangianScalarFieldSource()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::LagrangianSubScalarField>
Foam::totalNumberLagrangianScalarFieldSource::value
(
    const LagrangianInjection& injection,
    const LagrangianSubMesh& subMesh
) const
{
    // Calculate the necessary sizes
    tmp<LagrangianSubScalarField> size, v, m;
    calcSizes
    (
        injection, subMesh,
        size,
        haveVolume_, v,
        !haveVolume_, m
    );

    // Return the numbers that equalise the sizes on all parcels and recover
    // the specified total
    if (haveVolume_)
    {
        const dimensionedScalar V(dimVolume, volumeOrMass_);

        return V/size()/sum(v()/size());
    }
    else
    {
        const dimensionedScalar M(dimMass, volumeOrMass_);

        return M/size()/sum(m()/size());
    }
}


void Foam::totalNumberLagrangianScalarFieldSource::write(Ostream& os) const
{
    uniformSizeNumberLagrangianScalarFieldSource::write(os);

    if (haveVolume_)
    {
        writeEntry(os, "volume", volumeOrMass_);
    }
    else
    {
        writeEntry(os, "mass", volumeOrMass_);
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makeLagrangianTypeFieldSource
    (
        LagrangianScalarFieldSource,
        totalNumberLagrangianScalarFieldSource
    );
}

// ************************************************************************* //
