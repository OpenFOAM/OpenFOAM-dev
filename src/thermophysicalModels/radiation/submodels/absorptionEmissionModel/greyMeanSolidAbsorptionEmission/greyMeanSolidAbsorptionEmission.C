/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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

#include "greyMeanSolidAbsorptionEmission.H"
#include "addToRunTimeSelectionTable.H"
#include "unitConversion.H"
#include "extrapolatedCalculatedFvPatchFields.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace radiation
    {
        defineTypeNameAndDebug(greyMeanSolidAbsorptionEmission, 0);

        addToRunTimeSelectionTable
        (
            absorptionEmissionModel,
            greyMeanSolidAbsorptionEmission,
            dictionary
        );
    }
}

// * * * * * * * * * * * * * * Private Member Functions * * * * * * * * * * //

Foam::tmp<Foam::scalarField> Foam::radiation::
greyMeanSolidAbsorptionEmission::X(const word specie) const
{
    const volScalarField& T = thermo_.T();
    const volScalarField& p = thermo_.p();

    tmp<scalarField> tXj(new scalarField(T.primitiveField().size(), 0.0));
    scalarField& Xj = tXj.ref();

    tmp<scalarField> tRhoInv(new scalarField(T.primitiveField().size(), 0.0));
    scalarField& rhoInv = tRhoInv.ref();

    forAll(mixture_.Y(), specieI)
    {
        const scalarField& Yi = mixture_.Y()[specieI];

        forAll(rhoInv, iCell)
        {
            rhoInv[iCell] +=
                Yi[iCell]/mixture_.rho(specieI, p[iCell], T[iCell]);
        }
    }
    const scalarField& Yj = mixture_.Y(specie);
    const label mySpecieI = mixture_.species()[specie];
    forAll(Xj, iCell)
    {
        Xj[iCell] = Yj[iCell]/mixture_.rho(mySpecieI, p[iCell], T[iCell]);
    }

    return (Xj/rhoInv);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::radiation::greyMeanSolidAbsorptionEmission::
greyMeanSolidAbsorptionEmission
(
    const dictionary& dict,
    const fvMesh& mesh
)
:
    absorptionEmissionModel(dict, mesh),
    coeffsDict_((dict.optionalSubDict(typeName + "Coeffs"))),
    thermo_(mesh.lookupObject<solidThermo>(basicThermo::dictName)),
    speciesNames_(0),
    mixture_(dynamic_cast<const basicSpecieMixture&>(thermo_)),
    solidData_(mixture_.Y().size())
{
    if (!isA<basicSpecieMixture>(thermo_))
    {
        FatalErrorInFunction
            << "Model requires a multi-component thermo package"
            << abort(FatalError);
    }

    label nFunc = 0;
    const dictionary& functionDicts = dict.optionalSubDict(typeName + "Coeffs");

    forAllConstIter(dictionary, functionDicts, iter)
    {
        // safety:
        if (!iter().isDict())
        {
            continue;
        }
        const word& key = iter().keyword();
        if (!mixture_.contains(key))
        {
            WarningInFunction
                << " specie: " << key << " is not found in the solid mixture"
                << nl
                << " specie is the mixture are:" << mixture_.species() << nl
                << nl << endl;
        }
        speciesNames_.insert(key, nFunc);
        const dictionary& dict = iter().dict();
        dict.lookup("absorptivity") >> solidData_[nFunc][absorptivity];
        dict.lookup("emissivity") >> solidData_[nFunc][emissivity];

        nFunc++;
    }
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::radiation::greyMeanSolidAbsorptionEmission::
~greyMeanSolidAbsorptionEmission()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::radiation::greyMeanSolidAbsorptionEmission::
calc(const label propertyId) const
{
    tmp<volScalarField> ta
    (
        new volScalarField
        (
            IOobject
            (
                "a",
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh(),
            dimensionedScalar("a", dimless/dimLength, 0.0),
            extrapolatedCalculatedFvPatchVectorField::typeName
        )
    );

    scalarField& a = ta.ref().primitiveFieldRef();

    forAllConstIter(HashTable<label>, speciesNames_, iter)
    {
        if (mixture_.contains(iter.key()))
        {
            a += solidData_[iter()][propertyId]*X(iter.key());
        }
    }

    ta.ref().correctBoundaryConditions();
    return ta;
}


Foam::tmp<Foam::volScalarField>
Foam::radiation::greyMeanSolidAbsorptionEmission::eCont
(
    const label bandI
) const
{
   return calc(emissivity);
}


Foam::tmp<Foam::volScalarField>
Foam::radiation::greyMeanSolidAbsorptionEmission::aCont
(
    const label bandI
) const
{
   return calc(absorptivity);
}

// ************************************************************************* //
