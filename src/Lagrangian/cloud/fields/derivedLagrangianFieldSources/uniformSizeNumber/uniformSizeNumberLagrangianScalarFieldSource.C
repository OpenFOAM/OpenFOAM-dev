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

#include "uniformSizeNumberLagrangianScalarFieldSource.H"
#include "massive.H"

// * * * * * * * * * * * * * Static Member Data  * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(uniformSizeNumberLagrangianScalarFieldSource, 0);
}


namespace Foam
{
    template<>
    const char* NamedEnum
    <
        uniformSizeNumberLagrangianScalarFieldSource::uniformSize,
        4
    >::names[] = {"number", "surfaceArea", "volume", "mass"};
}


const Foam::NamedEnum
<
    Foam::uniformSizeNumberLagrangianScalarFieldSource::uniformSize,
    4
> Foam::uniformSizeNumberLagrangianScalarFieldSource::uniformSizeNames_;


// * * * * * * * * * * * *  Protected Member Functions * * * * * * * * * * * //

void Foam::uniformSizeNumberLagrangianScalarFieldSource::calcSizes
(
    const LagrangianInjection& injection,
    const LagrangianSubMesh& subMesh,
    tmp<LagrangianSubScalarField>& size,
    const bool needv,
    tmp<LagrangianSubScalarField>& v,
    const bool needm,
    tmp<LagrangianSubScalarField>& m
) const
{
    const bool actuallyNeedv = needv || uniformSize_ == uniformSize::volume;
    const bool actuallyNeedm = needm || uniformSize_ == uniformSize::mass;

    const Foam::clouds::shaped& shapedCloud =
        this->cloud<Foam::clouds::shaped>(injection, subMesh);

    // Evaluate the volume and mass (as necessary) of the created particles
    if (actuallyNeedv)
    {
        v = shapedCloud.v(injection, subMesh);
    }
    if (actuallyNeedm)
    {
        const clouds::massive& massiveCloud =
            this->cloud<clouds::massive>(injection, subMesh);

        m = massiveCloud.m(injection, subMesh);
    }

    // Create a field proportional to the desired uniform size
    switch (uniformSize_)
    {
        case uniformSize::number:
        {
            size =
                LagrangianSubScalarField::New
                (
                    "1",
                    subMesh,
                    dimensionedScalar(dimless, 1)
                );
            break;
        }
        case uniformSize::surfaceArea:
        {
            size = shapedCloud.a(injection, subMesh);
            break;
        }
        case uniformSize::volume:
        {
            size = v();
            break;
        }
        case uniformSize::mass:
        {
            size = m();
            break;
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::uniformSizeNumberLagrangianScalarFieldSource::
uniformSizeNumberLagrangianScalarFieldSource
(
    const regIOobject& iIo,
    const dictionary& dict
)
:
    LagrangianScalarFieldSource(iIo, dict),
    CloudLagrangianFieldSource<scalar>(*this),
    uniformSize_(uniformSizeNames_[dict.lookup<word>("uniformSize")])
{}


Foam::uniformSizeNumberLagrangianScalarFieldSource::
uniformSizeNumberLagrangianScalarFieldSource
(
    const uniformSizeNumberLagrangianScalarFieldSource& field,
    const regIOobject& iIo
)
:
    LagrangianScalarFieldSource(field, iIo),
    CloudLagrangianFieldSource<scalar>(*this),
    uniformSize_(field.uniformSize_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::uniformSizeNumberLagrangianScalarFieldSource::
~uniformSizeNumberLagrangianScalarFieldSource()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::label Foam::uniformSizeNumberLagrangianScalarFieldSource::sampleQ() const
{
    switch (uniformSize_)
    {
        case uniformSize::number:
            return 0;
        case uniformSize::surfaceArea:
            return 2;
        case uniformSize::volume:
            return 3;
        case uniformSize::mass:
            return 3;
    }

    return -labelMax;
}


void Foam::uniformSizeNumberLagrangianScalarFieldSource::write
(
    Ostream& os
) const
{
    LagrangianScalarFieldSource::write(os);

    writeEntry(os, "uniformSize", uniformSizeNames_[uniformSize_]);
}


// ************************************************************************* //
