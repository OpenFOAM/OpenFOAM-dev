/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2017-2018 OpenFOAM Foundation
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

#include "constantSurfaceTension.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace surfaceTensionModels
{
    defineTypeNameAndDebug(constant, 0);
    addToRunTimeSelectionTable(surfaceTensionModel, constant, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::surfaceTensionModels::constant::constant
(
    const dictionary& dict,
    const fvMesh& mesh
)
:
    surfaceTensionModel(mesh),
    sigma_("sigma", dimensionSet(1, 0, -2, 0, 0), dict)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::surfaceTensionModels::constant::~constant()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::surfaceTensionModels::constant::sigma() const
{
    return volScalarField::New(sigma_.name(), mesh_, sigma_);
}


bool Foam::surfaceTensionModels::constant::readDict(const dictionary& dict)
{
    // Handle sub-dictionary format as a special case
    if (dict.isDict("sigma"))
    {
        dict.subDict("sigma").lookup("sigma") >> sigma_;
    }
    else
    {
        dict.lookup("sigma") >> sigma_;
    }

    return true;
}


bool Foam::surfaceTensionModels::constant::writeData(Ostream& os) const
{
    if (surfaceTensionModel::writeData(os))
    {
        os  << sigma_ << token::END_STATEMENT << nl;
        return os.good();
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
