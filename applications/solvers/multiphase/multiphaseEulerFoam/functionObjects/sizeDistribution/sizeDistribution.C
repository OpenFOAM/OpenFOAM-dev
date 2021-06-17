/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2017-2021 OpenFOAM Foundation
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

#include "sizeDistribution.H"
#include "Time.H"
#include "fvMesh.H"
#include "addToRunTimeSelectionTable.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(sizeDistribution, 0);
    addToRunTimeSelectionTable(functionObject, sizeDistribution, dictionary);
}
}

template<>
const char*
Foam::NamedEnum
<
    Foam::functionObjects::sizeDistribution::functionType,
    4
>::names[] =
{
    "moments",
    "standardDeviation",
    "number",
    "volume"
};

const Foam::NamedEnum
<
    Foam::functionObjects::sizeDistribution::functionType,
    4
> Foam::functionObjects::sizeDistribution::functionTypeNames_;

template<>
const char*
Foam::NamedEnum
<
    Foam::functionObjects::sizeDistribution::coordinateType,
    4
>::names[] =
{
    "volume",
    "area",
    "diameter",
    "projectedAreaDiameter"
};

const Foam::NamedEnum
<
    Foam::functionObjects::sizeDistribution::coordinateType,
    4
> Foam::functionObjects::sizeDistribution::coordinateTypeNames_;

using Foam::constant::mathematical::pi;


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

Foam::tmp<Foam::scalarField>
Foam::functionObjects::sizeDistribution::filterField
(
    const scalarField& field
) const
{
    if (isNull(cellIDs()))
    {
        return field;
    }
    else
    {
        return tmp<scalarField>(new scalarField(field, cellIDs()));
    }
}


void Foam::functionObjects::sizeDistribution::correctVolAverages()
{
    forAll(N_, i)
    {
        const Foam::diameterModels::sizeGroup& fi = popBal_.sizeGroups()[i];

        scalarField Ni(filterField(fi*fi.phase()/fi.x()));
        scalarField V(filterField(mesh_.V()));
        scalarField ai(filterField(fi.a()));
        scalarField di(Ni.size());

        switch (coordinateType_)
        {
            case ctProjectedAreaDiameter:
            {
                di = sqrt(ai/pi);

                break;
            }
            default:
            {
                di = fi.d();

                break;
            }
        }

        N_[i] = gSum(V*Ni)/this->V();
        a_[i] = gSum(V*ai)/this->V();
        d_[i] = gSum(V*di)/this->V();
    }

    forAll(bins_, i)
    {
        const Foam::diameterModels::sizeGroup& fi = popBal_.sizeGroups()[i];

        bins_[i] = point(fi.x().value(), a_[i], d_[i]);
    }
}


void Foam::functionObjects::sizeDistribution::writeMoments()
{
    logFiles::write();

    Log << "    writing moments of size distribution." << endl;

    if (Pstream::master())
    {
        writeTime(file());
    }

    for (label k = 0; k <= maxOrder_; k++)
    {
        scalar result = 0;

        forAll(N_, i)
        {
            result += pow(bins_[i][binCmpt_], k)*N_[i];
        }

        if (Pstream::master())
        {
            file() << tab << result;
        }
    }

    if (Pstream::master())
    {
        file() << endl;
    }
}


void Foam::functionObjects::sizeDistribution::writeStdDev()
{
    logFiles::write();

    Log << "    writing standard deviation of size distribution."
        << endl;

    if (Pstream::master())
    {
        writeTime(file());
    }

    scalar stdDev = 0;
    scalar mean = 0;
    scalar var = 0;

    if(sum(N_) != 0)
    {
        if (geometric_)
        {
            mean = exp(sum(Foam::log(bins_.component(binCmpt_))*N_/sum(N_)));

            var =
                sum(sqr(Foam::log(bins_.component(binCmpt_)) - Foam::log(mean))
               *N_/sum(N_));

            stdDev = exp(sqrt(var));
        }
        else
        {
            mean = sum(bins_.component(binCmpt_)*N_/sum(N_));

            var = sum(sqr(bins_.component(binCmpt_) - mean)*N_/sum(N_));

            stdDev = sqrt(var);
        }
    }

    if (Pstream::master())
    {
        file() << tab << stdDev << tab << mean << tab << var << endl;
    }
}


void Foam::functionObjects::sizeDistribution::writeDistribution()
{
    scalarField result(N_);

    switch (functionType_)
    {
        case ftNumber:
        {
            Log << "    writing number distribution. "
                << endl;

            break;
        }
        case ftVolume:
        {
            Log << "    writing volume distribution. "
                << endl;

            result *= bins_.component(0);

            break;
        }
        default:
        {
            break;
        }

    }

    if (normalise_)
    {
        if(sum(result) != 0)
        {
            result /= sum(result);
        }

    }

    if (densityFunction_)
    {
        List<scalar> bndrs(N_.size() + 1);

        bndrs.first() = bins_.first()[binCmpt_];
        bndrs.last() = bins_.last()[binCmpt_];

        for (label i = 1; i < N_.size(); i++)
        {
            bndrs[i] = (bins_[i][binCmpt_] + bins_[i-1][binCmpt_])/2.0;
        }

        forAll(result, i)
        {
            if (geometric_)
            {
                result[i] /=
                    (Foam::log(bndrs[i+1]) - Foam::log(bndrs[i]));
            }
            else
            {
                result[i] /= (bndrs[i+1] - bndrs[i]);
            }
        }
    }

    if (Pstream::master())
    {
        const coordSet coords
        (
            "sizeDistribution",
            "xyz",
            bins_,
            mag(bins_)
        );

        writeGraph(coords, functionTypeNames_[functionType_], result);
    }
}


void Foam::functionObjects::sizeDistribution::writeFileHeader
(
    const label i
)
{
    volRegion::writeFileHeader(*this, file());

    writeHeaderValue
    (
        file(),
        "Coordinate",
        word(coordinateTypeNames_[coordinateType_])
    );

    word str("Time");

    switch (functionType_)
    {
        case ftMoments:
        {
            for (label k = 0; k <= maxOrder_; k++)
            {
                str += (" k=" + std::to_string(k));
            }

            break;
        }

        case ftStdDev:
        {
            str += " standardDeviation mean variance";

            break;
        }

        default:
        {
            break;
        }
    }

    writeCommented(file(), str);

    file() << endl;
}


void Foam::functionObjects::sizeDistribution::writeGraph
(
    const coordSet& coords,
    const word& functionTypeName,
    const scalarField& values
)
{
    const wordList functionTypeNames(1, functionTypeName);

    fileName outputPath = file_.baseTimeDir();

    mkDir(outputPath);
    OFstream graphFile
    (
        outputPath/(this->name() + ".dat")
    );

    volRegion::writeFileHeader(file_, graphFile);

    file_.writeCommented(graphFile, "Volume area diameter " + functionTypeName);

    if (densityFunction_)
    {
        graphFile << "Density";
    }
    else
    {
        graphFile << "Concentration";
    }

    graphFile << endl;

    List<const scalarField*> yPtrs(1);
    yPtrs[0] = &values;
    scalarFormatter_().write(coords, functionTypeNames, yPtrs, graphFile);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::sizeDistribution::sizeDistribution
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    volRegion(fvMeshFunctionObject::mesh_, dict),
    logFiles(obr_, name),
    mesh_(fvMeshFunctionObject::mesh_),
    file_(obr_, name),
    popBal_
    (
        obr_.lookupObject<Foam::diameterModels::populationBalanceModel>
        (
            dict.lookup("populationBalance")
        )
    ),
    functionType_(functionTypeNames_.read(dict.lookup("functionType"))),
    coordinateType_(coordinateTypeNames_.read(dict.lookup("coordinateType"))),
    N_(popBal_.sizeGroups().size(), 0),
    a_(popBal_.sizeGroups().size(), 0),
    d_(popBal_.sizeGroups().size(), 0),
    bins_(N_.size()),
    binCmpt_(0)
{
    read(dict);

    switch (coordinateType_)
    {
        case ctVolume:
        {
            binCmpt_ = 0;

            break;
        }
        case ctArea:
        {
            binCmpt_ = 1;

            break;
        }
        case ctDiameter:
        {
            binCmpt_ = 2;

            break;
        }
        case ctProjectedAreaDiameter:
        {
            binCmpt_ = 2;

            break;
        }
    }

    scalarFormatter_ = setWriter<scalar>::New("raw");
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::sizeDistribution::~sizeDistribution()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::sizeDistribution::read(const dictionary& dict)
{
    fvMeshFunctionObject::read(dict);

    normalise_ = dict.lookupOrDefault<Switch>("normalise", false);
    densityFunction_ = dict.lookupOrDefault<Switch>("densityFunction", false);
    geometric_ = dict.lookupOrDefault<Switch>("geometric", false);
    maxOrder_ = dict.lookupOrDefault("maxOrder", 3);

    resetName(name());

    return false;
}


bool Foam::functionObjects::sizeDistribution::execute()
{
    return true;
}


bool Foam::functionObjects::sizeDistribution::end()
{
    return true;
}


bool Foam::functionObjects::sizeDistribution::write()
{
    Log << type() << " " << name() << " write:" << nl;

    correctVolAverages();

    switch (functionType_)
    {
        case ftMoments:
        {
            writeMoments();

            break;
        }

        case ftStdDev:
        {
            writeStdDev();

            break;
        }

        default:
        {
            writeDistribution();

            break;
        }
    }

    Log << endl;

    return true;
}


// ************************************************************************* //
