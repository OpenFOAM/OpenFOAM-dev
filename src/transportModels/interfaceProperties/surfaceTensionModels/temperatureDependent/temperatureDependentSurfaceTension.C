/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2017 OpenFOAM Foundation
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

#include "temperatureDependentSurfaceTension.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace surfaceTensionModels
{
    defineTypeNameAndDebug(temperatureDependent, 0);
    addToRunTimeSelectionTable
    (
        surfaceTensionModel,
        temperatureDependent,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::surfaceTensionModels::temperatureDependent::temperatureDependent
(
    const dictionary& dict,
    const fvMesh& mesh
)
:
    surfaceTensionModel(mesh),
    TName_(dict.lookupOrDefault<word>("T", "T")),
    sigma_(Function1<scalar>::New("sigma", dict))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::surfaceTensionModels::temperatureDependent::~temperatureDependent()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::surfaceTensionModels::temperatureDependent::sigma() const
{
    tmp<volScalarField> tsigma
    (
        new volScalarField
        (
            IOobject
            (
                "sigma",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh_,
            dimSigma
        )
    );
    volScalarField& sigma = tsigma.ref();

    const volScalarField& T = mesh_.lookupObject<volScalarField>(TName_);

    sigma.field() = sigma_->value(T.field());

    volScalarField::Boundary& sigmaBf = sigma.boundaryFieldRef();
    const volScalarField::Boundary& TBf = T.boundaryField();

    forAll(sigmaBf, patchi)
    {
        sigmaBf[patchi] = sigma_->value(TBf[patchi]);
    }

    return tsigma;
}


bool Foam::surfaceTensionModels::temperatureDependent::readDict
(
    const dictionary& dict
)
{
    const dictionary& sigmaDict = surfaceTensionModel::sigmaDict(dict);

    TName_ = sigmaDict.lookupOrDefault<word>("T", "T");
    sigma_ = Function1<scalar>::New("sigma", sigmaDict);

    return true;
}


bool Foam::surfaceTensionModels::temperatureDependent::writeData
(
    Ostream& os
) const
{
    if (surfaceTensionModel::writeData(os))
    {
        os  << sigma_() << token::END_STATEMENT << nl;
        return os.good();
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
