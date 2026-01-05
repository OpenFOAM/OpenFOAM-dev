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

#include "basicLagrangianThermo.H"
#include "calculatedLagrangianPatchFields.H"
#include "densityLagrangianScalarFieldSource.H"
#include "specificHeatCapacityLagrangianScalarFieldSource.H"
#include "thermalConductivityLagrangianScalarFieldSource.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(basicLagrangianThermo, 0);
    defineRunTimeSelectionTable(basicLagrangianThermo, LagrangianMesh);
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

Foam::wordList Foam::basicLagrangianThermo::eBoundaryTypes() const
{
    const LagrangianScalarDynamicField::Boundary& Tbf = T().boundaryField();

    wordList eBt = Tbf.types();

    // !!! There are no boundary conditions for temperature in Lagrangian other
    // than calculated and constraint types. So the types are the same for the
    // temperature and energy fields. If in future boundary conditions are
    // added that affect the temperature (e.g., some sort of wall-heat transfer
    // model) then corresponding energy conditions will also be needed and this
    // function will need to translate between the two.

    return eBt;
}


Foam::wordList Foam::basicLagrangianThermo::eBoundaryBaseTypes() const
{
    const LagrangianScalarDynamicField::Boundary& Tbf = T().boundaryField();

    wordList eBbt(Tbf.size(), word::null);

    // !!! There is no "overrides constraint" mechanism in Lagrangian at
    // present. So there is currently nothing to be done here.

    return eBbt;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::basicLagrangianThermo::implementation::implementation
(
    const dictionary& dict,
    const LagrangianMesh& mesh,
    const word& phaseName
)
:
    mesh_(mesh),
    phaseName_(phaseName),
    T_
    (
        IOobject
        (
            IOobject::groupName("T", phaseName),
            mesh.time().name(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ),
    rho_
    (
        IOobject
        (
            IOobject::groupName("rho", phaseName),
            mesh.time().name(),
            mesh
        ),
        mesh,
        dimensionedScalar("NaN", dimDensity, NaN),
        wordList
        (
            mesh.boundary().size(),
            calculatedLagrangianPatchScalarField::typeName
        ),
        wordList::null(),
        sourcesTypes<densityLagrangianScalarFieldSource>(T_),
        T_.sources().errorLocation()
    ),
    Cv_
    (
        IOobject
        (
            IOobject::groupName("Cv", phaseName),
            mesh.time().name(),
            mesh
        ),
        mesh,
        dimensionedScalar("NaN", dimSpecificHeatCapacity, NaN),
        wordList
        (
            mesh.boundary().size(),
            calculatedLagrangianPatchScalarField::typeName
        ),
        wordList::null(),
        sourcesTypes<specificHeatCapacityLagrangianScalarFieldSource>(T_),
        T_.sources().errorLocation()
    ),
    kappa_
    (
        IOobject
        (
            IOobject::groupName("kappa", phaseName),
            mesh.time().name(),
            mesh
        ),
        mesh,
        dimensionedScalar("NaN", dimThermalConductivity, NaN),
        wordList
        (
            mesh.boundary().size(),
            calculatedLagrangianPatchScalarField::typeName
        ),
        wordList::null(),
        sourcesTypes<thermalConductivityLagrangianScalarFieldSource>(T_),
        T_.sources().errorLocation()
    )
{}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::basicLagrangianThermo>
Foam::basicLagrangianThermo::New
(
    const LagrangianMesh& mesh,
    const word& phaseName
)
{
    return New<basicLagrangianThermo>(mesh, phaseName);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::basicLagrangianThermo::~basicLagrangianThermo()
{}


Foam::basicLagrangianThermo::implementation::~implementation()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

const Foam::LagrangianMesh&
Foam::basicLagrangianThermo::implementation::mesh() const
{
    return mesh_;
}


const Foam::word&
Foam::basicLagrangianThermo::implementation::phaseName() const
{
    return phaseName_;
}


const Foam::LagrangianScalarDynamicField&
Foam::basicLagrangianThermo::implementation::T() const
{
    return T_;
}


Foam::LagrangianScalarDynamicField&
Foam::basicLagrangianThermo::implementation::T()
{
    return T_;
}


const Foam::LagrangianScalarDynamicField&
Foam::basicLagrangianThermo::implementation::rho() const
{
    return rho_;
}


Foam::LagrangianScalarDynamicField&
Foam::basicLagrangianThermo::implementation::rho()
{
    return rho_;
}


const Foam::LagrangianScalarDynamicField&
Foam::basicLagrangianThermo::implementation::Cv() const
{
    return Cv_;
}


const Foam::LagrangianScalarDynamicField&
Foam::basicLagrangianThermo::implementation::kappa() const
{
    return kappa_;
}


void Foam::basicLagrangianThermo::implementation::read(const dictionary&)
{}


// ************************************************************************* //
