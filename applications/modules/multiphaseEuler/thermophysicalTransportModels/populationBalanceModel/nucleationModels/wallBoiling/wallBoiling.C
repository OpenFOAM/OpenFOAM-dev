/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2018-2024 OpenFOAM Foundation
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

#include "wallBoiling.H"
#include "addToRunTimeSelectionTable.H"
#include "alphatWallBoilingWallFunctionFvPatchScalarField.H"
#include "phaseSystem.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace diameterModels
{
namespace nucleationModels
{
    defineTypeNameAndDebug(wallBoiling, 0);
    addToRunTimeSelectionTable
    (
        nucleationModel,
        wallBoiling,
        dictionary
    );
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::diameterModels::nucleationModels::wallBoiling::
wallBoiling
(
    const populationBalanceModel& popBal,
    const dictionary& dict
)
:
    nucleationModel(popBal, dict)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::diameterModels::nucleationModels::wallBoiling::precompute()
{
    const volScalarField& alphat =
        popBal_.mesh().lookupObject<volScalarField>
        (
            IOobject::groupName("alphat", popBal_.continuousPhase().name())
        );

    const volScalarField::Boundary& alphatBf = alphat.boundaryField();

    typedef compressible::alphatWallBoilingWallFunctionFvPatchScalarField
        alphatWallBoilingWallFunction;

    forAll(alphatBf, patchi)
    {
        if (isA<alphatWallBoilingWallFunction>(alphatBf[patchi]))
        {
            const alphatWallBoilingWallFunction& alphatw =
                refCast<const alphatWallBoilingWallFunction>(alphatBf[patchi]);

            const scalarField& dDep = alphatw.dDeparture();

            if (min(dDep) < velGroup_.sizeGroups().first().dSph().value())
            {
                Warning
                    << "Minimum departure diameter " << min(dDep)
                    << " m outside of range ["
                    << velGroup_.sizeGroups().first().dSph().value() << ", "
                    << velGroup_.sizeGroups().last().dSph().value() << "] m"
                    << " at patch " << alphatw.patch().name()
                    << endl
                    << "    The nucleation rate in populationBalance "
                    << popBal_.name() << " is set to zero." << endl
                    << "    Adjust discretisation over property space to"
                    << " suppress this warning."
                    << endl;
            }
            else if (max(dDep) > velGroup_.sizeGroups().last().dSph().value())
            {
                Warning
                    << "Maximum departure diameter " << max(dDep)
                    << " m outside of range ["
                    << velGroup_.sizeGroups().first().dSph().value() << ", "
                    << velGroup_.sizeGroups().last().dSph().value() << "] m"
                    << " at patch " << alphatw.patch().name()
                    << endl
                    << "    The nucleation rate in populationBalance "
                    << popBal_.name() << " is set to zero." << endl
                    << "    Adjust discretisation over property space to"
                    << " suppress this warning."
                    << endl;
            }
        }
    }
}


void Foam::diameterModels::nucleationModels::wallBoiling::addToNucleationRate
(
    volScalarField& nucleationRate,
    const label i
)
{
    const sizeGroup& fi = popBal_.sizeGroups()[i];
    const phaseModel& phase = fi.phase();
    const volScalarField& rho = phase.rho();

    const volScalarField& alphat =
        popBal_.mesh().lookupObject<volScalarField>
        (
            IOobject::groupName("alphat", popBal_.continuousPhase().name())
        );

    const volScalarField::Boundary& alphatBf = alphat.boundaryField();

    typedef compressible::alphatWallBoilingWallFunctionFvPatchScalarField
        alphatWallBoilingWallFunction;

    forAll(alphatBf, patchi)
    {
        if (!isA<alphatWallBoilingWallFunction>(alphatBf[patchi])) continue;

        const alphatWallBoilingWallFunction& alphatw =
            refCast<const alphatWallBoilingWallFunction>(alphatBf[patchi]);

        const scalarField& dmdt = alphatw.dmdtf();
        const scalarField& dDep = alphatw.dDeparture();

        const labelList& faceCells = alphatw.patch().faceCells();

        dimensionedScalar unitLength("unitLength", dimLength, 1);

        forAll(alphatw, facei)
        {
            if (dmdt[facei] > small)
            {
                const label faceCelli = faceCells[facei];

                nucleationRate[faceCelli] +=
                    popBal_.eta
                    (
                        i,
                        fi.x()/pow3(fi.dSph())*pow3(dDep[facei]*unitLength)
                    ).value()
                   *dmdt[facei]/rho[faceCelli]/fi.x().value();
            }
        }
    }
}


// ************************************************************************* //
