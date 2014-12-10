/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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

#include "basicThermo.H"
#include "zeroGradientFvPatchFields.H"
#include "fixedEnergyFvPatchScalarField.H"
#include "gradientEnergyFvPatchScalarField.H"
#include "mixedEnergyFvPatchScalarField.H"
#include "fixedJumpFvPatchFields.H"
#include "fixedJumpAMIFvPatchFields.H"
#include "energyJumpFvPatchScalarField.H"
#include "energyJumpAMIFvPatchScalarField.H"


/* * * * * * * * * * * * * * * private static data * * * * * * * * * * * * * */

namespace Foam
{
    defineTypeNameAndDebug(basicThermo, 0);
    defineRunTimeSelectionTable(basicThermo, fvMesh);
}

const Foam::word Foam::basicThermo::dictName("thermophysicalProperties");


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

Foam::wordList Foam::basicThermo::heBoundaryBaseTypes()
{
    const volScalarField::GeometricBoundaryField& tbf =
        this->T_.boundaryField();

    wordList hbt(tbf.size(), word::null);

    forAll(tbf, patchi)
    {
        if (isA<fixedJumpFvPatchScalarField>(tbf[patchi]))
        {
            const fixedJumpFvPatchScalarField& pf =
                dynamic_cast<const fixedJumpFvPatchScalarField&>(tbf[patchi]);

            hbt[patchi] = pf.interfaceFieldType();
        }
        else if (isA<fixedJumpAMIFvPatchScalarField>(tbf[patchi]))
        {
            const fixedJumpAMIFvPatchScalarField& pf =
                dynamic_cast<const fixedJumpAMIFvPatchScalarField&>
                (
                    tbf[patchi]
                );

            hbt[patchi] = pf.interfaceFieldType();
        }
    }

    return hbt;
}


Foam::wordList Foam::basicThermo::heBoundaryTypes()
{
    const volScalarField::GeometricBoundaryField& tbf =
        this->T_.boundaryField();

    wordList hbt = tbf.types();

    forAll(tbf, patchi)
    {
        if (isA<fixedValueFvPatchScalarField>(tbf[patchi]))
        {
            hbt[patchi] = fixedEnergyFvPatchScalarField::typeName;
        }
        else if
        (
            isA<zeroGradientFvPatchScalarField>(tbf[patchi])
         || isA<fixedGradientFvPatchScalarField>(tbf[patchi])
        )
        {
            hbt[patchi] = gradientEnergyFvPatchScalarField::typeName;
        }
        else if (isA<mixedFvPatchScalarField>(tbf[patchi]))
        {
            hbt[patchi] = mixedEnergyFvPatchScalarField::typeName;
        }
        else if (isA<fixedJumpFvPatchScalarField>(tbf[patchi]))
        {
            hbt[patchi] = energyJumpFvPatchScalarField::typeName;
        }
        else if (isA<fixedJumpAMIFvPatchScalarField>(tbf[patchi]))
        {
            hbt[patchi] = energyJumpAMIFvPatchScalarField::typeName;
        }
        else if (tbf[patchi].type() == "energyRegionCoupledFvPatchScalarField")
        {
            hbt[patchi] = "energyRegionCoupledFvPatchScalarField";
        }
    }

    return hbt;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::volScalarField& Foam::basicThermo::lookupOrConstruct
(
    const fvMesh& mesh,
    const char* name
) const
{
    if (!mesh.objectRegistry::foundObject<volScalarField>(name))
    {
        volScalarField* fPtr
        (
            new volScalarField
            (
                IOobject
                (
                    name,
                    mesh.time().timeName(),
                    mesh,
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh
            )
        );

        // Transfer ownership of this object to the objectRegistry
        fPtr->store(fPtr);
    }

    return const_cast<volScalarField&>
    (
        mesh.objectRegistry::lookupObject<volScalarField>(name)
    );
}


Foam::basicThermo::basicThermo
(
    const fvMesh& mesh,
    const word& phaseName
)
:
    IOdictionary
    (
        IOobject
        (
            phasePropertyName(dictName, phaseName),
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    ),

    phaseName_(phaseName),

    p_(lookupOrConstruct(mesh, "p")),

    T_
    (
        IOobject
        (
            phasePropertyName("T"),
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ),

    alpha_
    (
        IOobject
        (
            phasePropertyName("thermo:alpha"),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionSet(1, -1, -1, 0, 0)
    ),

    dpdt_(lookupOrDefault<Switch>("dpdt", true))
{}


Foam::basicThermo::basicThermo
(
    const fvMesh& mesh,
    const dictionary& dict,
    const word& phaseName
)
:
    IOdictionary
    (
        IOobject
        (
            phasePropertyName(dictName, phaseName),
            mesh.time().constant(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        dict
    ),

    phaseName_(phaseName),

    p_(lookupOrConstruct(mesh, "p")),

    T_
    (
        IOobject
        (
            phasePropertyName("T"),
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ),

    alpha_
    (
        IOobject
        (
            phasePropertyName("thermo:alpha"),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionSet(1, -1, -1, 0, 0)
    )
{}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::basicThermo> Foam::basicThermo::New
(
    const fvMesh& mesh,
    const word& phaseName
)
{
    return New<basicThermo>(mesh, phaseName);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::basicThermo::~basicThermo()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::basicThermo& Foam::basicThermo::lookupThermo
(
    const fvPatchScalarField& pf
)
{
    if (pf.db().foundObject<basicThermo>(dictName))
    {
        return pf.db().lookupObject<basicThermo>(dictName);
    }
    else
    {
        HashTable<const basicThermo*> thermos =
            pf.db().lookupClass<basicThermo>();

        for
        (
            HashTable<const basicThermo*>::iterator iter = thermos.begin();
            iter != thermos.end();
            ++iter
        )
        {
            if
            (
                &(iter()->he().dimensionedInternalField())
              == &(pf.dimensionedInternalField())
            )
            {
                return *iter();
            }
        }
    }

    return pf.db().lookupObject<basicThermo>(dictName);
}


void Foam::basicThermo::validate
(
    const string& app,
    const word& a
) const
{
    if (!(he().name() == phasePropertyName(a)))
    {
        FatalErrorIn(app)
            << "Supported energy type is " << phasePropertyName(a)
            << ", thermodynamics package provides " << he().name()
            << exit(FatalError);
    }
}

void Foam::basicThermo::validate
(
    const string& app,
    const word& a,
    const word& b
) const
{
    if
    (
       !(
            he().name() == phasePropertyName(a)
         || he().name() == phasePropertyName(b)
        )
    )
    {
        FatalErrorIn(app)
            << "Supported energy types are " << phasePropertyName(a)
            << " and " << phasePropertyName(b)
            << ", thermodynamics package provides " << he().name()
            << exit(FatalError);
    }
}

void Foam::basicThermo::validate
(
    const string& app,
    const word& a,
    const word& b,
    const word& c
) const
{
    if
    (
       !(
            he().name() == phasePropertyName(a)
         || he().name() == phasePropertyName(b)
         || he().name() == phasePropertyName(c)
        )
    )
    {
        FatalErrorIn(app)
            << "Supported energy types are " << phasePropertyName(a)
            << ", " << phasePropertyName(b)
            << " and " << phasePropertyName(c)
            << ", thermodynamics package provides " << he().name()
            << exit(FatalError);
    }
}

void Foam::basicThermo::validate
(
    const string& app,
    const word& a,
    const word& b,
    const word& c,
    const word& d
) const
{
    if
    (
       !(
            he().name() == phasePropertyName(a)
         || he().name() == phasePropertyName(b)
         || he().name() == phasePropertyName(c)
         || he().name() == phasePropertyName(d)
        )
    )
    {
        FatalErrorIn(app)
            << "Supported energy types are " << phasePropertyName(a)
            << ", " << phasePropertyName(b)
            << ", " << phasePropertyName(c)
            << " and " << phasePropertyName(d)
            << ", thermodynamics package provides " << he().name()
            << exit(FatalError);
    }
}


Foam::wordList Foam::basicThermo::splitThermoName
(
    const word& thermoName,
    const int nCmpt
)
{
    wordList cmpts(nCmpt);

    string::size_type beg=0, end=0, endb=0, endc=0;
    int i = 0;

    while
    (
        (endb = thermoName.find('<', beg)) != string::npos
     || (endc = thermoName.find(',', beg)) != string::npos
    )
    {
        if (endb == string::npos)
        {
            end = endc;
        }
        else if ((endc = thermoName.find(',', beg)) != string::npos)
        {
            end = min(endb, endc);
        }
        else
        {
            end = endb;
        }

        if (beg < end)
        {
            cmpts[i] = thermoName.substr(beg, end-beg);
            cmpts[i++].replaceAll(">","");
        }
        beg = end + 1;
    }

    if (beg < thermoName.size())
    {
        cmpts[i] = thermoName.substr(beg, string::npos);
        cmpts[i++].replaceAll(">","");
    }

    return cmpts;
}


Foam::volScalarField& Foam::basicThermo::p()
{
    return p_;
}


const Foam::volScalarField& Foam::basicThermo::p() const
{
    return p_;
}


const Foam::volScalarField& Foam::basicThermo::T() const
{
    return T_;
}


Foam::volScalarField& Foam::basicThermo::T()
{
    return T_;
}


const Foam::volScalarField& Foam::basicThermo::alpha() const
{
    return alpha_;
}


const Foam::scalarField& Foam::basicThermo::alpha(const label patchi) const
{
    return alpha_.boundaryField()[patchi];
}


bool Foam::basicThermo::read()
{
    return regIOobject::read();
}


// ************************************************************************* //
