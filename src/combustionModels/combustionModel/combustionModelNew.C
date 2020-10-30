/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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

#include "combustionModel.H"
#include "noCombustion.H"

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::combustionModel> Foam::combustionModel::New
(
    const fluidReactionThermo& thermo,
    const compressibleMomentumTransportModel& turb,
    const word& combustionProperties
)
{
    IOobject combIO
    (
        IOobject
        (
            thermo.phasePropertyName(combustionProperties),
            thermo.T().mesh().time().constant(),
            thermo.T().mesh(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        )
    );

    word modelType(combustionModels::noCombustion::typeName);
    if (combIO.typeHeaderOk<IOdictionary>(false))
    {
        IOdictionary(combIO).lookup(combustionModel::typeName) >> modelType;
    }
    else
    {
        Info<< "Combustion model not active: "
            << thermo.phasePropertyName(combustionProperties)
            << " not found" << endl;
    }

    Info<< "Selecting combustion model " << modelType << endl;

    const wordList cmpts2(basicThermo::splitThermoName(modelType, 2));
    const wordList cmpts3(basicThermo::splitThermoName(modelType, 3));
    if (cmpts2.size() == 2 || cmpts3.size() == 3)
    {
        modelType = cmpts2.size() ? cmpts2[0] : cmpts3[0];

        WarningInFunction
            << "Template parameters are no longer required when selecting a "
            << combustionModel::typeName << ". This information is now "
            << "obtained directly from the thermodynamics. Actually selecting "
            << "combustion model " << modelType << "." << endl;
    }

    typename dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(modelType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown " << combustionModel::typeName << " type "
            << modelType << nl << nl
            << "Valid " << combustionModel::typeName << " types are:" << nl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);

        const wordList names(dictionaryConstructorTablePtr_->sortedToc());
    }

    return autoPtr<combustionModel>
    (
        cstrIter()(modelType, thermo, turb, combustionProperties)
    );
}


// ************************************************************************* //
