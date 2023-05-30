/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2017-2023 OpenFOAM Foundation
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
#include "volFields.H"
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
    return
        mesh_.lookupObject<liquidThermo>
        (
             IOobject::groupName(physicalProperties::typeName, phaseName_)
        ).sigma();
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
