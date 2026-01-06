/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2026 OpenFOAM Foundation
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

#include "fluidLagrangianThermo.H"
#include "calculatedLagrangianPatchFields.H"
#include "pressureLagrangianScalarFieldSource.H"
#include "compressibilityLagrangianScalarFieldSource.H"
#include "dynamicViscosityLagrangianScalarFieldSource.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(fluidLagrangianThermo, 0);
    defineRunTimeSelectionTable(fluidLagrangianThermo, LagrangianMesh);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fluidLagrangianThermo::implementation::implementation
(
    const dictionary& dict,
    const LagrangianMesh& mesh,
    const word& phaseName,
    const LagrangianScalarDynamicField& T
)
:
    p_
    (
        IOobject
        (
            IOobject::groupName("p", phaseName),
            mesh.time().name(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("NaN", dimPressure, NaN),
        wordList
        (
            mesh.boundary().size(),
            calculatedLagrangianPatchScalarField::typeName
        ),
        wordList::null(),
        sourcesTypes<pressureLagrangianScalarFieldSource>(T),
        T.sources().errorLocation()
    ),
    pSourcePtr_
    (
        LagrangianScalarFieldSource::New
        (
            p_,
            dict.subDict("pressure")
        )
    ),
    psi_
    (
        IOobject
        (
            IOobject::groupName("psi", phaseName),
            mesh.time().name(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("NaN", dimDensity/dimPressure, NaN),
        wordList
        (
            mesh.boundary().size(),
            calculatedLagrangianPatchScalarField::typeName
        ),
        wordList::null(),
        sourcesTypes<compressibilityLagrangianScalarFieldSource>(T),
        T.sources().errorLocation()
    ),
    mu_
    (
        IOobject
        (
            IOobject::groupName("mu", phaseName),
            mesh.time().name(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("NaN", dimDynamicViscosity, NaN),
        wordList
        (
            mesh.boundary().size(),
            calculatedLagrangianPatchScalarField::typeName
        ),
        wordList::null(),
        sourcesTypes<dynamicViscosityLagrangianScalarFieldSource>(T),
        T.sources().errorLocation()
    )
{}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::fluidLagrangianThermo>
Foam::fluidLagrangianThermo::New
(
    const LagrangianMesh& mesh,
    const word& phaseName
)
{
    return basicLagrangianThermo::New<fluidLagrangianThermo>(mesh, phaseName);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fluidLagrangianThermo::~fluidLagrangianThermo()
{}


Foam::fluidLagrangianThermo::implementation::~implementation()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::fluidLagrangianThermo::implementation::initialise()
{
    // Initialise the pressure if it was not read in
    if (!p_.headerOk()) correctPressure(mesh().subAll());

    correct(mesh().subAll());
}


void Foam::fluidLagrangianThermo::implementation::correctPressure
(
    const LagrangianSubMesh& subMesh
)
{
    if (debug) InfoInFunction << endl;

    SubField<scalar> p = subMesh.sub(this->p_.primitiveFieldRef());

    // Use the source model to update the pressure on this sub-mesh
    p = pSourcePtr_->value(subMesh)().primitiveField();

    if (debug) Info<< "    Finished" << endl;
}


const Foam::LagrangianScalarDynamicField&
Foam::fluidLagrangianThermo::implementation::p() const
{
    return p_;
}


Foam::LagrangianScalarDynamicField&
Foam::fluidLagrangianThermo::implementation::p()
{
    return p_;
}


Foam::tmp<Foam::LagrangianSubScalarSubField>
Foam::fluidLagrangianThermo::implementation::p
(
    const LagrangianSubMesh& subMesh
) const
{
    return subMesh.sub(p_);
}


Foam::tmp<Foam::LagrangianSubScalarField>
Foam::fluidLagrangianThermo::implementation::p
(
    const LagrangianInjection& injection,
    const LagrangianSubMesh& subMesh
) const
{
    return pSourcePtr_->value(injection, subMesh);
}


const Foam::LagrangianScalarDynamicField&
Foam::fluidLagrangianThermo::implementation::psi() const
{
    return psi_;
}


const Foam::LagrangianScalarDynamicField&
Foam::fluidLagrangianThermo::implementation::mu() const
{
    return mu_;
}


// ************************************************************************* //
