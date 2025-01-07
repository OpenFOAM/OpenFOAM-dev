/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2025 OpenFOAM Foundation
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

#include "LagrangianDistribution.H"
#include "LagrangianFields.H"
#include "OSspecific.H"
#include "writeFile.H"
#include "unintegrable.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(LagrangianDistribution, 0);
    addToRunTimeSelectionTable
    (
        functionObject,
        LagrangianDistribution,
        dictionary
    );
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::functionObjects::LagrangianDistribution::readCoeffs
(
    const dictionary& dict
)
{
    // Read the fields
    const bool haveFields = dict.found("fields");
    const bool haveField = dict.found("field");
    if (haveFields == haveField)
    {
        FatalIOErrorInFunction(dict)
            << "keywords fields and field both "
            << (haveFields ? "" : "un") << "defined in "
            << "dictionary " << dict.name()
            << exit(FatalIOError);
    }
    else if (haveFields)
    {
        dict.lookup("fields") >> fields_;
    }
    else if (haveField)
    {
        fields_.resize(1);
        dict.lookup("field") >> fields_.first();
    }

    // Read the weight fields
    const bool haveWeightFields = dict.found("weightFields");
    const bool haveWeightField = dict.found("weightField");
    if (haveWeightFields && haveWeightField)
    {
        FatalIOErrorInFunction(dict)
            << "keywords weightFields and weightField both "
            << "defined in dictionary " << dict.name()
            << exit(FatalIOError);
    }
    else if (haveWeightFields)
    {
        dict.lookup("weightFields") >> weightFields_;
    }
    else if (haveWeightField)
    {
        weightFields_.resize(1);
        dict.lookup("weightField") >> weightFields_.first();
    }
    else
    {
        // No weights
        weightFields_.clear();
    }

    // The number of bins in the plot
    dict.lookup("nBins") >> nBins_;

    // Construct the formatter
    formatter_ = setWriter::New(dict.lookup("setFormat"), dict);
}


template<template<class> class GeoField>
bool Foam::functionObjects::LagrangianDistribution::multiplyWeight
(
    const word& weightFieldName,
    scalarField& weight
) const
{
    if (!mesh().foundObject<GeoField<scalar>>(weightFieldName)) return false;

    const GeoField<scalar>& w =
        mesh().lookupObject<GeoField<scalar>>(weightFieldName);

    weight *= w;

    return true;
}


void Foam::functionObjects::LagrangianDistribution::writeDistribution
(
    const scalarField& weight,
    const word& fieldName,
    const scalarField& field
)
{
    // Get the limits of the distribution
    const scalar x0 = gMin(field), xN = gMax(field);

    // Construct the limits of the bins
    scalarField x(nBins_ + 1);
    for (label nodei = 0; nodei <= nBins_; ++ nodei)
    {
        x[nodei] = x0 + nodei/scalar(nBins_)*(xN - x0);
    }

    // Populate the bins
    scalarField PDF(nBins_ + 1, scalar(0));
    forAll(field, i)
    {
        const scalar x = field[i];
        const scalar f = (x - x0)/max(xN - x0, rootVSmall);
        const scalar bini = min(floor(f*nBins_), nBins_ - 1);
        const scalar g = f*nBins_ - scalar(bini);
        PDF[bini] += weight[i]*(1 - g);
        PDF[bini + 1] += weight[i]*g;
    }

    // Synchronise
    Pstream::listCombineGather(PDF, plusEqOp<scalar>());
    Pstream::listCombineScatter(PDF);

    // Normalise and correct the ends, as they have half as many samples as the
    // interior points
    PDF /= sum(PDF)*(xN - x0)/nBins_;
    PDF.first() *= 2;
    PDF.last() *= 2;

    // Write
    if (Pstream::master())
    {
        const fileName outputPath =
            time_.globalPath()
           /writeFile::outputPrefix
           /(
                mesh().mesh().name() != polyMesh::defaultRegion
              ? mesh().mesh().name()
              : word::null
            )
           /name()
           /time_.name();

        mkDir(outputPath);

        formatter_->write
        (
            outputPath,
            fieldName,
            coordSet(true, fieldName, x),
            "PDF", PDF,
            "CDF", distributions::unintegrable::integrate(x, PDF)()
        );
    }
}


template<template<class> class GeoField, class Type>
bool Foam::functionObjects::LagrangianDistribution::writeDistribution
(
    const scalarField& weight,
    const word& fieldName
)
{
    if (!mesh().foundObject<GeoField<Type>>(fieldName)) return false;

    const typename GeoField<Type>::FieldType& field =
        mesh().lookupObject<GeoField<Type>>(fieldName).primitiveField();

    for (direction d = 0; d < pTraits<Type>::nComponents; ++ d)
    {
        writeDistribution
        (
            weight,
            fieldName
          + (word(pTraits<Type>::componentNames[d]).empty() ? "" : "_")
          + word(pTraits<Type>::componentNames[d]),
            field.component(d)()
        );
    }

    return true;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::LagrangianDistribution::LagrangianDistribution
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    LagrangianMeshFunctionObject(name, runTime, dict),
    fields_(),
    weightFields_(),
    nBins_(-1),
    formatter_(nullptr)
{
    readCoeffs(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::LagrangianDistribution::~LagrangianDistribution()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::LagrangianDistribution::read(const dictionary& dict)
{
    if (LagrangianMeshFunctionObject::read(dict))
    {
        readCoeffs(dict);
        return true;
    }
    else
    {
        return false;
    }
}


Foam::wordList Foam::functionObjects::LagrangianDistribution::fields() const
{
    wordList result(fields_);
    result.append(weightFields_);
    return result;
}


bool Foam::functionObjects::LagrangianDistribution::execute()
{
    return true;
}


bool Foam::functionObjects::LagrangianDistribution::write()
{
    // We can't construct a distribution without any samples, so quit if the
    // mesh is empty.
    if (returnReduce(mesh().size(), sumOp<label>()) == 0)
    {
        return true;
    }

    // Construct the weights
    scalarField weight(mesh().size(), 1);
    forAll(weightFields_, weightFieldi)
    {
        const word& weightFieldName = weightFields_[weightFieldi];

        if
        (
            !multiplyWeight<LagrangianField>(weightFieldName, weight)
         && !multiplyWeight<LagrangianDynamicField>(weightFieldName, weight)
         && !multiplyWeight<LagrangianInternalField>(weightFieldName, weight)
        )
        {
            FatalErrorInFunction
                << "Weight field " << weightFieldName << " was not found"
                << exit(FatalError);
        }
    }

    // Write the field values
    forAll(fields_, fieldi)
    {
        const word& fieldName = fields_[fieldi];

        #define WRITE_FIELD_VALUE(Type, GeoField) \
            && !writeDistribution<GeoField, Type>(weight, fieldName)

        if
        (
            true
            FOR_ALL_FIELD_TYPES(WRITE_FIELD_VALUE, LagrangianField)
            FOR_ALL_FIELD_TYPES(WRITE_FIELD_VALUE, LagrangianDynamicField)
            FOR_ALL_FIELD_TYPES(WRITE_FIELD_VALUE, LagrangianInternalField)
        )
        {
            cannotFindObject(fields_[fieldi]);
        }

        #undef WRITE_COLUMN_VALUE
    }

    return true;
}


// ************************************************************************* //
