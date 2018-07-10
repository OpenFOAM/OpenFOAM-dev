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

#include "wideBandAbsorptionEmission.H"
#include "addToRunTimeSelectionTable.H"
#include "basicSpecieMixture.H"
#include "unitConversion.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace radiation
    {
        defineTypeNameAndDebug(wideBandAbsorptionEmission, 0);

        addToRunTimeSelectionTable
        (
            absorptionEmissionModel,
            wideBandAbsorptionEmission,
            dictionary
        );
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::radiation::wideBandAbsorptionEmission::wideBandAbsorptionEmission
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
    Yj_(nSpecies_),
    totalWaveLength_(0)
{
    label nBand = 0;
    const dictionary& functionDicts = dict.optionalSubDict(typeName +"Coeffs");
    forAllConstIter(dictionary, functionDicts, iter)
    {
        // safety:
        if (!iter().isDict())
        {
            continue;
        }

        const dictionary& dict = iter().dict();
        dict.lookup("bandLimits") >> iBands_[nBand];
        dict.lookup("EhrrCoeff") >> iEhrrCoeffs_[nBand];
        totalWaveLength_ += iBands_[nBand][1] - iBands_[nBand][0];

        label nSpec = 0;
        const dictionary& specDicts = dict.subDict("species");
        forAllConstIter(dictionary, specDicts, iter)
        {
            const word& key = iter().keyword();
            if (nBand == 0)
            {
                speciesNames_.insert(key, nSpec);
            }
            else
            {
                if (!speciesNames_.found(key))
                {
                    FatalErrorInFunction
                        << "specie: " << key << " is not in all the bands"
                        << nl << exit(FatalError);
                }
            }
            coeffs_[nBand][nSpec].initialise(specDicts.subDict(key));
            nSpec++;
        }
        nBand++;
    }
    nBands_ = nBand;

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
                const label index =
                    lookUpTablePtr_().findFieldIndex(iter.key());

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
                << " there is no lookup table and the specie" << nl
                << iter.key() << nl
                << " is not found " << nl
                << exit(FatalError);

        }
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::radiation::wideBandAbsorptionEmission::~wideBandAbsorptionEmission()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::radiation::wideBandAbsorptionEmission::aCont(const label bandi) const
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
                "a",
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh(),
            dimensionedScalar("a", dimless/dimLength, 0)
        )
    );

    scalarField& a = ta.ref().primitiveFieldRef();

    forAll(a, celli)
    {
        forAllConstIter(HashTable<label>, speciesNames_, iter)
        {
            const label n = iter();
            scalar Xipi = 0;
            if (specieIndex_[n] != 0)
            {
                const volScalarField& ft =
                    mesh_.lookupObject<volScalarField>("ft");

                const List<scalar>& Ynft = lookUpTablePtr_().lookUp(ft[celli]);

                // moles*pressure [atm]
                Xipi = Ynft[specieIndex_[n]]*paToAtm(p[celli]);
            }
            else
            {
                scalar invWt = 0;
                forAll(mixture.Y(), s)
                {
                    invWt += mixture.Y(s)[celli]/mixture.W(s);
                }

                const label index = mixture.species()[iter.key()];

                const scalar Xk =
                    mixture.Y(index)[celli]/(mixture.W(index)*invWt);

                Xipi = Xk*paToAtm(p[celli]);
            }

            scalar Ti = T[celli];

            const absorptionCoeffs::coeffArray& b =
                coeffs_[bandi][n].coeffs(T[celli]);

            if (coeffs_[bandi][n].invTemp())
            {
                Ti = 1.0/T[celli];
            }

            a[celli]+=
                Xipi
               *(
                    ((((b[5]*Ti + b[4])*Ti + b[3])*Ti + b[2])*Ti + b[1])*Ti
                  + b[0]
                );
        }
    }

    return ta;
}


Foam::tmp<Foam::volScalarField>
Foam::radiation::wideBandAbsorptionEmission::eCont(const label bandi) const
{
    return aCont(bandi);
}


Foam::tmp<Foam::volScalarField>
Foam::radiation::wideBandAbsorptionEmission::ECont(const label bandi) const
{
    tmp<volScalarField> E
    (
        new volScalarField
        (
            IOobject
            (
                "E",
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh(),
            dimensionedScalar("E", dimMass/dimLength/pow3(dimTime), 0)
        )
    );

    if (mesh().foundObject<volScalarField>("Qdot"))
    {
        const volScalarField& Qdot =
            mesh().lookupObject<volScalarField>("Qdot");

        if (Qdot.dimensions() == dimEnergy/dimTime)
        {
            E.ref().primitiveFieldRef() =
                iEhrrCoeffs_[bandi]
               *Qdot.primitiveField()
               *(iBands_[bandi][1] - iBands_[bandi][0])
               /totalWaveLength_
               /mesh_.V();
        }
        else if (Qdot.dimensions() == dimEnergy/dimTime/dimVolume)
        {
            E.ref().primitiveFieldRef() =
                iEhrrCoeffs_[bandi]
               *Qdot.primitiveField()
               *(iBands_[bandi][1] - iBands_[bandi][0])
               /totalWaveLength_;
        }
        else
        {
            WarningInFunction
                << "Incompatible dimensions for Qdot field" << endl;
        }
    }

    return E;
}


void Foam::radiation::wideBandAbsorptionEmission::correct
(
    volScalarField& a,
    PtrList<volScalarField>& aLambda
) const
{
    a = dimensionedScalar("zero", dimless/dimLength, 0);

    for (label j=0; j<nBands_; j++)
    {
        aLambda[j].primitiveFieldRef() = this->a(j);

        a.primitiveFieldRef() +=
            aLambda[j].primitiveField()
           *(iBands_[j][1] - iBands_[j][0])
           /totalWaveLength_;
    }

}


// ************************************************************************* //
