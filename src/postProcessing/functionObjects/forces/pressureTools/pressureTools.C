/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012-2016 OpenFOAM Foundation
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

#include "pressureTools.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(pressureTools, 0);
    addToRunTimeSelectionTable(functionObject, pressureTools, dictionary);
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::word Foam::functionObjects::pressureTools::pName() const
{
    word fieldName = pName_;

    if (calcTotal_)
    {
        fieldName = "total(" + fieldName + ")";
    }
    else
    {
        fieldName = "static(" + fieldName + ")";
    }

    if (calcCoeff_)
    {
        fieldName = fieldName + "_coeff";
    }

    return fieldName;
}


Foam::dimensionedScalar Foam::functionObjects::pressureTools::rhoScale
(
    const volScalarField& p
) const
{
    if (p.dimensions() == dimPressure)
    {
        return dimensionedScalar("1", dimless,  1.0);
    }
    else
    {
        return dimensionedScalar("rhoRef", dimDensity, rhoInf_);
    }
}


Foam::tmp<Foam::volScalarField> Foam::functionObjects::pressureTools::rho
(
    const volScalarField& p
) const
{
    if (p.dimensions() == dimPressure)
    {
        return p.mesh().lookupObject<volScalarField>(rhoName_);
    }
    else
    {
        return
            tmp<volScalarField>
            (
                new volScalarField
                (
                    IOobject
                    (
                        "rho",
                        p.mesh().time().timeName(),
                        p.mesh(),
                        IOobject::NO_READ,
                        IOobject::NO_WRITE
                    ),
                    p.mesh(),
                    dimensionedScalar("zero", dimDensity, rhoInf_)
                )
            );
    }
}


Foam::dimensionedScalar Foam::functionObjects::pressureTools::pRef() const
{
    dimensionedScalar value("pRef", dimPressure, 0.0);

    if (calcTotal_)
    {
        value.value() += pRef_;
    }

    return value;
}


Foam::tmp<Foam::volScalarField> Foam::functionObjects::pressureTools::pDyn
(
    const volScalarField& p
) const
{
    tmp<volScalarField> tpDyn
    (
        new volScalarField
        (
            IOobject
            (
                "pDyn",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar("zero", dimPressure, 0.0)
        )
    );

    if (calcTotal_)
    {
        const volVectorField& U = obr_.lookupObject<volVectorField>(UName_);

        tpDyn.ref() == rho(p)*0.5*magSqr(U);
    }

    return tpDyn;
}


Foam::tmp<Foam::volScalarField>
Foam::functionObjects::pressureTools::convertToCoeff
(
    const volScalarField& p
) const
{
    tmp<volScalarField> tCoeff(p);

    if (calcCoeff_)
    {
        tCoeff.ref() -= dimensionedScalar("pInf", dimPressure, pInf_);

        const dimensionedScalar p0("p0", dimPressure, SMALL);
        const dimensionedVector U("U", dimVelocity, UInf_);
        const dimensionedScalar rho("rho", dimDensity, rhoInf_);

        tCoeff.ref() /= 0.5*rho*magSqr(U) + p0;
    }

    return tCoeff;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::pressureTools::pressureTools
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    pName_("p"),
    UName_("U"),
    rhoName_("rho"),
    calcTotal_(false),
    pRef_(0.0),
    calcCoeff_(false),
    pInf_(0.0),
    UInf_(Zero),
    rhoInf_(0.0)
{
    read(dict);

    dimensionSet pDims(dimPressure);

    if (calcCoeff_)
    {
        pDims /= dimPressure;
    }

    volScalarField* pPtr
    (
        new volScalarField
        (
            IOobject
            (
                pName(),
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar("0", pDims, 0.0)
        )
    );

    mesh_.objectRegistry::store(pPtr);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::pressureTools::~pressureTools()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::pressureTools::read(const dictionary& dict)
{
    dict.readIfPresent("pName", pName_);
    dict.readIfPresent("UName", UName_);
    dict.readIfPresent("rhoName", rhoName_);

    if (rhoName_ == "rhoInf")
    {
        dict.lookup("rhoInf") >> rhoInf_;
    }

    dict.lookup("calcTotal") >> calcTotal_;
    if (calcTotal_)
    {
        dict.lookup("pRef") >> pRef_;
    }

    dict.lookup("calcCoeff") >> calcCoeff_;
    if (calcCoeff_)
    {
        dict.lookup("pInf") >> pInf_;
        dict.lookup("UInf") >> UInf_;
        dict.lookup("rhoInf") >> rhoInf_;

        scalar zeroCheck = 0.5*rhoInf_*magSqr(UInf_) + pInf_;

        if (mag(zeroCheck) < ROOTVSMALL)
        {
            WarningInFunction
                << type() << " " << name() << ": "
                << "Coefficient calculation requested, but reference "
                << "pressure level is zero.  Please check the supplied "
                << "values of pInf, UInf and rhoInf" << endl;
        }
    }

    return true;
}


bool Foam::functionObjects::pressureTools::execute(const bool postProcess)
{
    const volScalarField& p = obr_.lookupObject<volScalarField>(pName_);

    volScalarField& pResult = const_cast<volScalarField&>
    (
        obr_.lookupObject<volScalarField>(pName())
    );

    pResult == convertToCoeff(rhoScale(p)*p + pDyn(p) + pRef());

    return true;
}


bool Foam::functionObjects::pressureTools::write(const bool postProcess)
{
    const volScalarField& pResult =
        obr_.lookupObject<volScalarField>(pName());

    Info<< type() << " " << name() << " output:" << nl
        << "    writing field " << pResult.name() << nl
        << endl;

    pResult.write();

    return true;
}


// ************************************************************************* //
