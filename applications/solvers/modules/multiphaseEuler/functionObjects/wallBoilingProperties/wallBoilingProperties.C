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
    volScalarField dDeparture
    (
        volScalarField::New
        (
            IOobject::groupName("dDeparture", phase_.name()),
            mesh_,
            dimensionedScalar(dimLength, 0)
        )
    );
    volScalarField fDeparture
    (
        volScalarField::New
        (
            IOobject::groupName("fDeparture", phase_.name()),
            mesh_,
            dimensionedScalar(inv(dimTime), 0)
        )
    );
    volScalarField nucleationSiteDensity
    (
        volScalarField::New
        (
            IOobject::groupName("nucleationSiteDensity", phase_.name()),
            mesh_,
            dimensionedScalar(inv(dimArea), 0)
        )
    );
    volScalarField wetFraction
    (
        volScalarField::New
        (
            IOobject::groupName("wetFraction", phase_.name()),
            mesh_,
            dimensionedScalar(dimless, 0)
        )
    );
    volScalarField qQuenching
    (
        volScalarField::New
        (
            IOobject::groupName("qQuenching", phase_.name()),
            mesh_,
            dimensionedScalar(dimEnergy*inv(dimTime*dimArea), 0)
        )
    );
    volScalarField qEvaporative
    (
        volScalarField::New
        (
            IOobject::groupName("qEvaporative", phase_.name()),
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

                    dDeparture.boundaryFieldRef()[patchi] =
                        alphatw.dDeparture();
                    fDeparture.boundaryFieldRef()[patchi] =
                        alphatw.fDeparture();
                    nucleationSiteDensity.boundaryFieldRef()[patchi] =
                        alphatw.nucleationSiteDensity();
                    wetFraction.boundaryFieldRef()[patchi] =
                        alphatw.wetFraction();
                    qQuenching.boundaryFieldRef()[patchi] =
                        alphatw.qQuenching();
                    qEvaporative.boundaryFieldRef()[patchi] =
                        alphatw.qEvaporative();
            }
        }
    }

    dDeparture.write();
    fDeparture.write();
    nucleationSiteDensity.write();
    wetFraction.write();
    qQuenching.write();
    qEvaporative.write();

    return true;
}


// ************************************************************************* //
