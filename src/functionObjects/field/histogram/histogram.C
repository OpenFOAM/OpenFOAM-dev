/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2021 OpenFOAM Foundation
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

#include "histogram.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(histogram, 0);
    addToRunTimeSelectionTable(functionObject, histogram, dictionary);
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::functionObjects::histogram::writeGraph
(
    const coordSet& coords,
    const word& fieldName,
    const scalarField& values
) const
{
    const wordList fieldNames(1, fieldName);

    fileName outputPath = file_.baseTimeDir();
    mkDir(outputPath);
    OFstream graphFile
    (
        outputPath/formatterPtr_().getFileName(coords, fieldNames)
    );

    Log << "    Writing histogram of " << fieldName
        << " to " << graphFile.name() << endl;

    List<const scalarField*> yPtrs(1);
    yPtrs[0] = &values;
    formatterPtr_().write(coords, fieldNames, yPtrs, graphFile);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::histogram::histogram
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    file_(obr_, name)
{
    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::histogram::~histogram()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::histogram::read(const dictionary& dict)
{
    dict.lookup("field") >> fieldName_;
    dict.lookup("max") >> max_;
    min_ = dict.lookupOrDefault<scalar>("min", 0);
    dict.lookup("nBins") >> nBins_;

    word format(dict.lookup("setFormat"));
    formatterPtr_ = setWriter<scalar>::New(format);

    return true;
}


Foam::wordList Foam::functionObjects::histogram::fields() const
{
    return wordList(fieldName_);
}


bool Foam::functionObjects::histogram::execute()
{
    return true;
}


bool Foam::functionObjects::histogram::write()
{
    Log << type() << " " << name() << " write:" << nl;

    autoPtr<volScalarField> fieldPtr;
    if (obr_.foundObject<volScalarField>(fieldName_))
    {
        Log << "    Looking up field " << fieldName_ << endl;
    }
    else
    {
        Log << "    Reading field " << fieldName_ << endl;
        fieldPtr.reset
        (
            new volScalarField
            (
                IOobject
                (
                    fieldName_,
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                ),
                mesh_
            )
        );
    }

    const volScalarField& field =
    (
         fieldPtr.valid()
       ? fieldPtr()
       : obr_.lookupObject<volScalarField>(fieldName_)
    );

    // Calculate the mid-points of bins for the graph axis
    pointField xBin(nBins_);
    const scalar delta = (max_- min_)/nBins_;

    scalar x = min_ + 0.5*delta;
    forAll(xBin, i)
    {
        xBin[i] = point(x, 0, 0);
        x += delta;
    }

    scalarField volFrac(nBins_, 0);
    const scalarField& V = mesh_.V();

    forAll(field, celli)
    {
        const label bini = (field[celli] - min_)/delta;
        if (bini >= 0 && bini < nBins_)
        {
            volFrac[bini] += V[celli];
        }
    }

    Pstream::listCombineGather(volFrac, plusEqOp<scalar>());

    if (Pstream::master())
    {
        const scalar sumVol = sum(volFrac);

        if (sumVol > small)
        {
            volFrac /= sumVol;

            const coordSet coords
            (
                "Volume_Fraction",
                "x",
                xBin,
                mag(xBin)
            );

            writeGraph(coords, field.name(), volFrac);
        }
    }

    return true;
}


// ************************************************************************* //
