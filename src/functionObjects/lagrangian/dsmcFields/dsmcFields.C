/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2022 OpenFOAM Foundation
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

#include "dsmcFields.H"
#include "volFields.H"
#include "dictionary.H"
#include "dsmcCloud.H"
#include "constants.H"
#include "addToRunTimeSelectionTable.H"

using namespace Foam::constant;

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(dsmcFields, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        dsmcFields,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::dsmcFields::dsmcFields
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict)
{
    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::dsmcFields::~dsmcFields()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::dsmcFields::read(const dictionary& dict)
{
    return true;
}


Foam::wordList Foam::functionObjects::dsmcFields::fields() const
{
    return wordList
    {
        "rhoNMean",
        "rhoMMean",
        "momentumMean",
        "linearKEMean",
        "internalEMean",
        "iDofMean",
        "fDMean"
    };
}


bool Foam::functionObjects::dsmcFields::execute()
{
    return true;
}


bool Foam::functionObjects::dsmcFields::write()
{
    word rhoNMeanName = "rhoNMean";
    word rhoMMeanName = "rhoMMean";
    word momentumMeanName = "momentumMean";
    word linearKEMeanName = "linearKEMean";
    word internalEMeanName = "internalEMean";
    word iDofMeanName = "iDofMean";
    word fDMeanName = "fDMean";

    const volScalarField& rhoNMean = obr_.lookupObject<volScalarField>
    (
        rhoNMeanName
    );

    const volScalarField& rhoMMean = obr_.lookupObject<volScalarField>
    (
        rhoMMeanName
    );

    const volVectorField& momentumMean = obr_.lookupObject<volVectorField>
    (
        momentumMeanName
    );

    const volScalarField& linearKEMean = obr_.lookupObject<volScalarField>
    (
        linearKEMeanName
    );

    const volScalarField& internalEMean = obr_.lookupObject<volScalarField>
    (
        internalEMeanName
    );

    const volScalarField& iDofMean = obr_.lookupObject<volScalarField>
    (
        iDofMeanName
    );

    const volVectorField& fDMean = obr_.lookupObject<volVectorField>
    (
        fDMeanName
    );

    if (min(mag(rhoNMean)).value() > vSmall)
    {
        Info<< "Calculating dsmcFields." << endl;

        Info<< "    Calculating UMean field." << endl;
        volVectorField UMean
        (
            IOobject
            (
                "UMean",
                obr_.time().name(),
                obr_,
                IOobject::NO_READ
            ),
            momentumMean/rhoMMean
        );

        Info<< "    Calculating translationalT field." << endl;
        volScalarField translationalT
        (
            IOobject
            (
                "translationalT",
                obr_.time().name(),
                obr_,
                IOobject::NO_READ
            ),

            2.0/(3.0*physicoChemical::k.value()*rhoNMean)
           *(linearKEMean - 0.5*rhoMMean*(UMean & UMean))
        );

        Info<< "    Calculating internalT field." << endl;
        volScalarField internalT
        (
            IOobject
            (
                "internalT",
                obr_.time().name(),
                obr_,
                IOobject::NO_READ
            ),
            (2.0/physicoChemical::k.value())*(internalEMean/iDofMean)
        );

        Info<< "    Calculating overallT field." << endl;
        volScalarField overallT
        (
            IOobject
            (
                "overallT",
                obr_.time().name(),
                obr_,
                IOobject::NO_READ
            ),
            2.0/(physicoChemical::k.value()*(3.0*rhoNMean + iDofMean))
           *(linearKEMean - 0.5*rhoMMean*(UMean & UMean) + internalEMean)
        );

        Info<< "    Calculating pressure field." << endl;
        volScalarField p
        (
            IOobject
            (
                "p",
                obr_.time().name(),
                obr_,
                IOobject::NO_READ
            ),
            physicoChemical::k.value()*rhoNMean*translationalT
        );

        volScalarField::Boundary& pBf = p.boundaryFieldRef();

        forAll(mesh_.boundaryMesh(), i)
        {
            const polyPatch& patch = mesh_.boundaryMesh()[i];

            if (isA<wallPolyPatch>(patch))
            {
                pBf[i] =
                    fDMean.boundaryField()[i]
                  & (patch.faceAreas()/patch.magFaceAreas());
            }
        }

        Info<< "    mag(UMean) max/min : "
            << max(mag(UMean)).value() << " "
            << min(mag(UMean)).value() << endl;

        Info<< "    translationalT max/min : "
            << max(translationalT).value() << " "
            << min(translationalT).value() << endl;

        Info<< "    internalT max/min : "
            << max(internalT).value() << " "
            << min(internalT).value() << endl;

        Info<< "    overallT max/min : "
            << max(overallT).value() << " "
            << min(overallT).value() << endl;

        Info<< "    p max/min : "
            << max(p).value() << " "
            << min(p).value() << endl;

        UMean.write();

        translationalT.write();

        internalT.write();

        overallT.write();

        p.write();

        Info<< "dsmcFields written." << nl << endl;

        return true;
    }
    else
    {
        Info<< "Small value (" << min(mag(rhoNMean))
            << ") found in rhoNMean field. "
            << "Not calculating dsmcFields to avoid division by zero."
            << endl;

        return false;
    }
}


// ************************************************************************* //
