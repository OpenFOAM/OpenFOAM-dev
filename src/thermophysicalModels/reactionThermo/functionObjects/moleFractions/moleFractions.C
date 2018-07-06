/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2018 OpenFOAM Foundation
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

#include "moleFractions.H"
#include "basicThermo.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class ThermoType>
void Foam::moleFractions<ThermoType>::calculateMoleFractions()
{
    const ThermoType& thermo =
        mesh_.lookupObject<ThermoType>(basicThermo::dictName);

    const PtrList<volScalarField>& Y = thermo.composition().Y();

    const volScalarField W(thermo.W());

    forAll(Y, i)
    {
        const dimensionedScalar Wi
        (
            "W",
            dimMass/dimMoles,
            thermo.composition().W(i)
        );

        X_[i] = W*Y[i]/Wi;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ThermoType>
Foam::moleFractions<ThermoType>::moleFractions
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict)
{
    if (mesh_.foundObject<ThermoType>(basicThermo::dictName))
    {
        const ThermoType& thermo =
            mesh_.lookupObject<ThermoType>(basicThermo::dictName);

        const PtrList<volScalarField>& Y = thermo.composition().Y();

        X_.setSize(Y.size());

        forAll(Y, i)
        {
            X_.set
            (
                i,
                new volScalarField
                (
                    IOobject
                    (
                        "X_" + Y[i].name(),
                        mesh_.time().timeName(),
                        mesh_,
                        IOobject::NO_READ,
                        IOobject::AUTO_WRITE
                    ),
                    mesh_,
                    dimensionedScalar("X", dimless, 0)
                )
            );
        }

        calculateMoleFractions();
    }
    else
    {
        FatalErrorInFunction
            << "Cannot find thermodynamics model of type "
            << ThermoType::typeName
            << exit(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class ThermoType>
Foam::moleFractions<ThermoType>::~moleFractions()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ThermoType>
bool Foam::moleFractions<ThermoType>::read
(
    const dictionary& dict
)
{
    return true;
}


template<class ThermoType>
bool Foam::moleFractions<ThermoType>::execute()
{
    calculateMoleFractions();
    return true;
}


template<class ThermoType>
bool Foam::moleFractions<ThermoType>::write()
{
    return true;
}


// ************************************************************************* //
