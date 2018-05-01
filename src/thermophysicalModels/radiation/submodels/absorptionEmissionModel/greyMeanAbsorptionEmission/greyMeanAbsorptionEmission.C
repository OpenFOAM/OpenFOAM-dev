/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
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

#include "greyMeanAbsorptionEmission.H"
#include "addToRunTimeSelectionTable.H"
#include "unitConversion.H"
#include "extrapolatedCalculatedFvPatchFields.H"
#include "basicSpecieMixture.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace radiation
    {
        defineTypeNameAndDebug(greyMeanAbsorptionEmission, 0);

        addToRunTimeSelectionTable
        (
            absorptionEmissionModel,
            greyMeanAbsorptionEmission,
            dictionary
        );
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::radiation::greyMeanAbsorptionEmission::greyMeanAbsorptionEmission
(
    const dictionary& dict,
    const fvMesh& mesh
)
:
    absorptionEmissionModel(dict, mesh),
    coeffsDict_((dict.optionalSubDict(typeName + "Coeffs"))),
    speciesNames_(0),
    specieIndex_(label(0)),
    lookUpTablePtr_(),
    thermo_(mesh.lookupObject<fluidThermo>(basicThermo::dictName)),
    EhrrCoeff_(readScalar(coeffsDict_.lookup("EhrrCoeff"))),
    Yj_(nSpecies_)
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
        speciesNames_.insert(key, nFunc);
        const dictionary& dict = iter().dict();
        coeffs_[nFunc].initialise(dict);
        nFunc++;
    }

    if (coeffsDict_.found("lookUpTableFileName"))
    {
        const word name = coeffsDict_.lookup("lookUpTableFileName");
        if (name != "none")
        {
            lookUpTablePtr_.set
            (
                new interpolationLookUpTable<scalar>
                (
                    fileName(coeffsDict_.lookup("lookUpTableFileName")),
                    mesh.time().constant(),
                    mesh
                )
            );

            if (!mesh.foundObject<volScalarField>("ft"))
            {
                FatalErrorInFunction
                    << "specie ft is not present to use with "
                    << "lookUpTableFileName " << nl
                    << exit(FatalError);
            }
        }
    }

    // Check that all the species on the dictionary are present in the
    // look-up table and save the corresponding indices of the look-up table

    label j = 0;
    forAllConstIter(HashTable<label>, speciesNames_, iter)
    {
        if (!lookUpTablePtr_.empty())
        {
            if (lookUpTablePtr_().found(iter.key()))
            {
                label index = lookUpTablePtr_().findFieldIndex(iter.key());

                Info<< "specie: " << iter.key() << " found on look-up table "
                    << " with index: " << index << endl;

                specieIndex_[iter()] = index;
            }
            else if (mesh.foundObject<volScalarField>(iter.key()))
            {
                Yj_.set(j, &mesh.lookupObjectRef<volScalarField>(iter.key()));
                specieIndex_[iter()] = 0;
                j++;
                Info<< "specie: " << iter.key() << " is being solved" << endl;
            }
            else
            {
                FatalErrorInFunction
                    << "specie: " << iter.key()
                    << " is neither in look-up table: "
                    << lookUpTablePtr_().tableName()
                    << " nor is being solved" << nl
                    << exit(FatalError);
            }
        }
        else if (mesh.foundObject<volScalarField>(iter.key()))
        {
            Yj_.set(j, &mesh.lookupObjectRef<volScalarField>(iter.key()));
            specieIndex_[iter()] = 0;
            j++;
        }
        else
        {
            FatalErrorInFunction
                << " there is not lookup table and the specie" << nl
                << iter.key() << nl
                << " is not found " << nl
                << exit(FatalError);

        }
    }
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::radiation::greyMeanAbsorptionEmission::~greyMeanAbsorptionEmission()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::radiation::greyMeanAbsorptionEmission::aCont(const label bandI) const
{
    const basicSpecieMixture& mixture =
        dynamic_cast<const basicSpecieMixture&>(thermo_);

    const volScalarField& T = thermo_.T();
    const volScalarField& p = thermo_.p();


    tmp<volScalarField> ta
    (
        new volScalarField
        (
            IOobject
            (
                "aCont" + name(bandI),
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

    forAll(a, celli)
    {
        forAllConstIter(HashTable<label>, speciesNames_, iter)
        {
            label n = iter();
            scalar Xipi = 0.0;
            if (specieIndex_[n] != 0)
            {
                // Specie found in the lookUpTable.
                const volScalarField& ft =
                    mesh_.lookupObject<volScalarField>("ft");

                const List<scalar>& Ynft = lookUpTablePtr_().lookUp(ft[celli]);
                // moles x pressure [atm]
                Xipi = Ynft[specieIndex_[n]]*paToAtm(p[celli]);
            }
            else
            {
                scalar invWt = 0.0;
                forAll(mixture.Y(), s)
                {
                    invWt += mixture.Y(s)[celli]/mixture.W(s);
                }

                label index = mixture.species()[iter.key()];
                scalar Xk = mixture.Y(index)[celli]/(mixture.W(index)*invWt);

                Xipi = Xk*paToAtm(p[celli]);
            }

            const absorptionCoeffs::coeffArray& b = coeffs_[n].coeffs(T[celli]);

            scalar Ti = T[celli];
            // negative temperature exponents
            if (coeffs_[n].invTemp())
            {
                Ti = 1.0/T[celli];
            }
            a[celli] +=
                Xipi
               *(
                    ((((b[5]*Ti + b[4])*Ti + b[3])*Ti + b[2])*Ti + b[1])*Ti
                  + b[0]
                );
        }
    }
    ta.ref().correctBoundaryConditions();
    return ta;
}


Foam::tmp<Foam::volScalarField>
Foam::radiation::greyMeanAbsorptionEmission::eCont(const label bandI) const
{
   return aCont(bandI);
}


Foam::tmp<Foam::volScalarField>
Foam::radiation::greyMeanAbsorptionEmission::ECont(const label bandI) const
{
    tmp<volScalarField> E
    (
        new volScalarField
        (
            IOobject
            (
                "ECont" + name(bandI),
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar("E", dimMass/dimLength/pow3(dimTime), 0.0)
        )
    );

    if (mesh_.foundObject<volScalarField>("Qdot"))
    {
        const volScalarField& Qdot =
            mesh_.lookupObject<volScalarField>("Qdot");

        if (Qdot.dimensions() == dimEnergy/dimTime)
        {
            E.ref().primitiveFieldRef() = EhrrCoeff_*Qdot/mesh_.V();
        }
        else if (Qdot.dimensions() == dimEnergy/dimTime/dimVolume)
        {
            E.ref().primitiveFieldRef() = EhrrCoeff_*Qdot;
        }
        else
        {
            if (debug)
            {
                WarningInFunction
                    << "Incompatible dimensions for Qdot field" << endl;
            }
        }
    }
    else
    {
        WarningInFunction
          << "Qdot field not found in mesh" << endl;
    }

    return E;
}


// ************************************************************************* //
