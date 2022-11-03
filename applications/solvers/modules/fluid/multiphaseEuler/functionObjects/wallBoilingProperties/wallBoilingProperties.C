/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2022 OpenFOAM Foundation
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

#include "wallBoilingProperties.H"
#include "addToRunTimeSelectionTable.H"
#include "alphatWallBoilingWallFunctionFvPatchScalarField.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(wallBoilingProperties, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        wallBoilingProperties,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::wallBoilingProperties::wallBoilingProperties
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    phase_
    (
        mesh_.lookupObject<phaseModel>
        (
            IOobject::groupName("alpha", dict.lookup("phase"))
        )
    ),
    fluid_(mesh_.lookupObject<phaseSystem>("phaseProperties"))
{
    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::wallBoilingProperties::~wallBoilingProperties()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::wallBoilingProperties::read(const dictionary& dict)
{
    fvMeshFunctionObject::read(dict);

    return true;
}


bool Foam::functionObjects::wallBoilingProperties::execute()
{
    return true;
}


bool Foam::functionObjects::wallBoilingProperties::write()
{
    volScalarField dDepartureField
    (
        volScalarField::New
        (
            IOobject::groupName("dDeparture", phase_.name()),
            mesh_,
            dimensionedScalar(dimLength, 0)
        )
    );
    volScalarField fDepartureField
    (
        volScalarField::New
        (
            IOobject::groupName("fDeparture", phase_.name()),
            mesh_,
            dimensionedScalar(inv(dimTime), 0)
        )
    );
    volScalarField nucSiteDensityField
    (
        volScalarField::New
        (
            IOobject::groupName("nucleationSiteDensity", phase_.name()),
            mesh_,
            dimensionedScalar(inv(dimArea), 0)
        )
    );
    volScalarField fLiquidField
    (
        volScalarField::New
        (
            IOobject::groupName("fLiquid", phase_.name()),
            mesh_,
            dimensionedScalar(dimless, 0)
        )
    );
    volScalarField quenchingHeatFluxField
    (
        volScalarField::New
        (
            IOobject::groupName("quenchingHeatFlux", phase_.name()),
            mesh_,
            dimensionedScalar(dimEnergy*inv(dimTime*dimArea), 0)
        )
    );
    volScalarField evaporativeHeatFluxField
    (
        volScalarField::New
        (
            IOobject::groupName("evaporativeHeatFlux", phase_.name()),
            mesh_,
            dimensionedScalar(dimEnergy*inv(dimTime*dimArea), 0)
        )
    );

    typedef compressible::alphatWallBoilingWallFunctionFvPatchScalarField
        alphatWallBoilingWallFunction;

    const word alphatName =
        IOobject::groupName("alphat", phase_.name());

    if (phase_.mesh().foundObject<volScalarField>(alphatName))
    {
        const volScalarField& alphat =
            phase_.mesh().lookupObject<volScalarField>(alphatName);

            const volScalarField::Boundary& alphatBf = alphat.boundaryField();

            forAll(alphatBf, patchi)
            {
                if (isA<alphatWallBoilingWallFunction>(alphatBf[patchi]))
                {
                    const alphatWallBoilingWallFunction& alphatw =
                        refCast
                        <
                            const alphatWallBoilingWallFunction
                        >(alphatBf[patchi]);

                    dDepartureField.boundaryFieldRef()[patchi] =
                        alphatw.dDeparture();
                    fDepartureField.boundaryFieldRef()[patchi] =
                        alphatw.depFrequency();
                    nucSiteDensityField.boundaryFieldRef()[patchi] =
                        alphatw.nucSiteDensity();
                    fLiquidField.boundaryFieldRef()[patchi] =
                        alphatw.wallLiquidFraction();
                    quenchingHeatFluxField.boundaryFieldRef()[patchi] =
                        alphatw.quenching();
                    evaporativeHeatFluxField.boundaryFieldRef()[patchi] =
                        alphatw.evaporative();
            }
        }
    }

    dDepartureField.write();
    fDepartureField.write();
    nucSiteDensityField.write();
    fLiquidField.write();
    quenchingHeatFluxField.write();
    evaporativeHeatFluxField.write();

    return true;
}


// ************************************************************************* //
