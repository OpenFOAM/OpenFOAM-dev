/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2017-2020 OpenFOAM Foundation
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

#include "liquidPropertiesSurfaceTension.H"
#include "liquidThermo.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace surfaceTensionModels
{
    defineTypeNameAndDebug(liquidProperties, 0);
    addToRunTimeSelectionTable
    (
        surfaceTensionModel,
        liquidProperties,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::surfaceTensionModels::liquidProperties::liquidProperties
(
    const dictionary& dict,
    const fvMesh& mesh
)
:
    surfaceTensionModel(mesh),
    phaseName_(dict.lookup("phase"))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::surfaceTensionModels::liquidProperties::~liquidProperties()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::surfaceTensionModels::liquidProperties::sigma() const
{
    const heRhoThermopureMixtureliquidProperties& thermo =
        mesh_.lookupObject<heRhoThermopureMixtureliquidProperties>
        (
             IOobject::groupName(basicThermo::dictName, phaseName_)
        );

    tmp<volScalarField> tsigma
    (
        volScalarField::New
        (
            "sigma",
            mesh_,
            dimSigma
        )
    );
    volScalarField& sigma = tsigma.ref();

    const volScalarField& T = thermo.T();
    const volScalarField& p = thermo.p();

    volScalarField::Internal& sigmai = sigma;
    const volScalarField::Internal& pi = p;
    const volScalarField::Internal& Ti = T;

    forAll(sigmai, celli)
    {
        const heRhoThermopureMixtureliquidProperties::mixtureType& liquid =
            thermo.cellMixture(celli);

        sigmai[celli] = liquid.properties().sigma(pi[celli], Ti[celli]);
    }

    volScalarField::Boundary& sigmaBf = sigma.boundaryFieldRef();
    const volScalarField::Boundary& pBf = p.boundaryField();
    const volScalarField::Boundary& TBf = T.boundaryField();

    forAll(sigmaBf, patchi)
    {
        scalarField& sigmaPf = sigmaBf[patchi];
        const scalarField& pPf = pBf[patchi];
        const scalarField& TPf = TBf[patchi];

        forAll(sigmaPf, facei)
        {
            const heRhoThermopureMixtureliquidProperties::mixtureType& liquid =
                thermo.patchFaceMixture(patchi, facei);

            sigmaPf[facei] = liquid.properties().sigma(pPf[facei], TPf[facei]);
        }
    }

    return tsigma;
}


bool Foam::surfaceTensionModels::liquidProperties::readDict
(
    const dictionary& dict
)
{
    return true;
}


bool Foam::surfaceTensionModels::liquidProperties::writeData
(
    Ostream& os
) const
{
    if (surfaceTensionModel::writeData(os))
    {
        return os.good();
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
