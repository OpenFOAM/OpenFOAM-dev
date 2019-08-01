/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2013-2019 OpenFOAM Foundation
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

#include "mixtureFraction.H"
#include "reactingMixture.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

template<class ThermoType>
const Foam::reactingMixture<ThermoType>&
Foam::radiationModels::sootModels::mixtureFraction<ThermoType>::checkThermo
(
    const fluidThermo& thermo
)
{
    if (isA<reactingMixture<ThermoType>>(thermo))
    {
        return dynamic_cast<const reactingMixture<ThermoType>& >
        (
            thermo
        );
    }
    else
    {
        FatalErrorInFunction
            << "Inconsistent thermo package for " << thermo.type()
            << "Please select a thermo package based on "
            << "reactingMixture" << exit(FatalError);

        return dynamic_cast<const reactingMixture<ThermoType>& >
        (
            thermo
        );
    }

}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ThermoType>
Foam::radiationModels::sootModels::mixtureFraction<ThermoType>::mixtureFraction
(
    const dictionary& dict,
    const fvMesh& mesh,
    const word& modelType
)
:
    sootModel(dict, mesh, modelType),
    soot_
    (
        IOobject
        (
            "soot",
            mesh_.time().timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),
    coeffsDict_(dict.subOrEmptyDict(modelType + "Coeffs")),
    nuSoot_(readScalar(coeffsDict_.lookup("nuSoot"))),
    Wsoot_(readScalar(coeffsDict_.lookup("Wsoot"))),
    sootMax_(-1),
    mappingFieldName_
    (
        coeffsDict_.lookupOrDefault<word>("mappingField", "none")
    ),
    mapFieldMax_(1),
    thermo_(mesh.lookupObject<fluidThermo>(basicThermo::dictName)),
    mixture_(checkThermo(thermo_))
{
    const Reaction<ThermoType>& reaction = mixture_.reactions()[0];

    scalar totalMol = 0;
    forAll(reaction.rhs(), i)
    {
        const scalar stoichCoeff = reaction.rhs()[i].stoichCoeff;
        totalMol += mag(stoichCoeff);
    }

    totalMol += nuSoot_;

    scalarList Xi(reaction.rhs().size());

    scalar Wm = 0;
    forAll(reaction.rhs(), i)
    {
        const label speciei = reaction.rhs()[i].index;
        const scalar stoichCoeff = reaction.rhs()[i].stoichCoeff;
        Xi[i] = mag(stoichCoeff)/totalMol;
        Wm += Xi[i]*mixture_.speciesData()[speciei].W();
    }

    scalarList Yprod0(mixture_.species().size(), 0.0);

    forAll(reaction.rhs(), i)
    {
        const label speciei = reaction.rhs()[i].index;
        Yprod0[speciei] = mixture_.speciesData()[speciei].W()/Wm*Xi[i];
    }

    const scalar XSoot = nuSoot_/totalMol;
    Wm += XSoot*Wsoot_;

    sootMax_ = XSoot*Wsoot_/Wm;

    Info << "Maximum soot mass concentrations: " << sootMax_ << nl;

    if (mappingFieldName_ == "none")
    {
        const label index = reaction.rhs()[0].index;
        mappingFieldName_ = mixture_.Y(index).name();
    }

    const label mapFieldIndex = mixture_.species()[mappingFieldName_];

    mapFieldMax_ = Yprod0[mapFieldIndex];
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class ThermoType>
Foam::radiationModels::sootModels::mixtureFraction<ThermoType>::
~mixtureFraction()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ThermoType>
void Foam::radiationModels::sootModels::mixtureFraction<ThermoType>::correct()
{
    const volScalarField& mapField =
        mesh_.lookupObject<volScalarField>(mappingFieldName_);

    soot_ = sootMax_*(mapField/mapFieldMax_);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
