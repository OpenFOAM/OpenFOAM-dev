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

#include "sizeDistribution.H"
#include "sizeGroup.H"
#include "addToRunTimeSelectionTable.H"

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
    Foam::functionObjects::sizeDistribution::regionTypes,
    2
>::names[] = {"cellZone", "all"};

template<>
const char*
Foam::NamedEnum
<
    Foam::functionObjects::sizeDistribution::functionTypes,
    4
>::names[] =
    {
        "numberDensity",
        "volumeDensity",
        "numberConcentration",
        "moments"
    };

template<>
const char*
Foam::NamedEnum
<
    Foam::functionObjects::sizeDistribution::abszissaTypes,
    2
>::names[] = {"diameter", "volume"};

const Foam::NamedEnum
<
    Foam::functionObjects::sizeDistribution::regionTypes,
    2
> Foam::functionObjects::sizeDistribution::regionTypeNames_;

const Foam::NamedEnum
<
    Foam::functionObjects::sizeDistribution::functionTypes,
    4
> Foam::functionObjects::sizeDistribution::functionTypeNames_;

const Foam::NamedEnum
<
    Foam::functionObjects::sizeDistribution::abszissaTypes,
    2
> Foam::functionObjects::sizeDistribution::abszissaTypeNames_;


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::functionObjects::sizeDistribution::initialise
(
    const dictionary& dict
)
{
    switch (functionType_)
    {
        case ftNdf:
        {
            break;
        }

        case ftVdf:
        {
            break;
        }

        case ftNc:
        {
            break;
        }

        case ftMom:
        {
            break;
        }

        default:
        {
            FatalErrorInFunction
               << "Unknown function type. Valid function types are:"
                << functionTypeNames_ << nl << exit(FatalError);
        }
    }

    switch (abszissaType_)
    {
        case atDiameter:
        {
            break;
        }

        case atVolume:
        {
            break;
        }

        default:
        {
            FatalErrorInFunction
               << "Unknown abszissa type. Valid abszissa types are:"
                << abszissaTypeNames_ << nl << exit(FatalError);
        }
    }

    setCellZoneCells();

    if (nCells_ == 0)
    {
        FatalErrorInFunction
            << type() << " " << name() << ": "
            << regionTypeNames_[regionType_] << "(" << regionName_ << "):" << nl
            << "    Region has no cells" << exit(FatalError);
    }

    volume_ = volume();

    Info<< type() << " " << name() << ":"
        << regionTypeNames_[regionType_] << "(" << regionName_ << "):" << nl
        << "    total cells  = " << nCells_ << nl
        << "    total volume = " << volume_
        << nl << endl;

    Info<< nl << endl;
}


void Foam::functionObjects::sizeDistribution::setCellZoneCells()
{
    switch (regionType_)
    {
        case rtCellZone:
        {
            dict().lookup("name") >> regionName_;

            label zoneId = mesh().cellZones().findZoneID(regionName_);

            if (zoneId < 0)
            {
                FatalErrorInFunction
                    << "Unknown cell zone name: " << regionName_
                    << ". Valid cell zones are: " << mesh().cellZones().names()
                    << nl << exit(FatalError);
            }

            cellId_ = mesh().cellZones()[zoneId];
            nCells_ = returnReduce(cellId_.size(), sumOp<label>());
            break;
        }

        case rtAll:
        {
            cellId_ = identity(mesh().nCells());
            nCells_ = returnReduce(cellId_.size(), sumOp<label>());
            break;
        }

        default:
        {
            FatalErrorInFunction
               << "Unknown region type. Valid region types are:"
                << regionTypeNames_ << nl << exit(FatalError);
        }
    }

    if (debug)
    {
        Pout<< "Selected region size = " << cellId_.size() << endl;
    }
}


Foam::scalar Foam::functionObjects::sizeDistribution::volume() const
{
    return gSum(filterField(mesh().V()));
}


void Foam::functionObjects::sizeDistribution::combineFields(scalarField& field)
{
    List<scalarField> allValues(Pstream::nProcs());

    allValues[Pstream::myProcNo()] = field;

    Pstream::gatherList(allValues);

    if (Pstream::master())
    {
        field =
            ListListOps::combine<scalarField>
            (
                allValues,
                accessOp<scalarField>()
            );
    }
}


Foam::tmp<Foam::scalarField>
Foam::functionObjects::sizeDistribution::filterField
(
    const scalarField& field
) const
{
    return tmp<scalarField>(new scalarField(field, cellId_));
}


void Foam::functionObjects::sizeDistribution::writeFileHeader
(
    const label i
)
{
    OFstream& file = this->file();

    switch (functionType_)
    {
        case ftNdf:
        {
            writeHeader(file, "Number density function");
            break;
        }

        case ftVdf:
        {
            writeHeader(file, "Volume density function");
            break;
        }

        case ftNc:
        {
            writeHeader(file, "Number concentration");
            break;
        }

        case ftMom:
        {
            writeHeader(file, "Moments");
            break;
        }
    }

    switch (abszissaType_)
    {
        case atVolume:
        {
            writeCommented(file, "Time/volume");
            break;
        }

        case atDiameter:
        {
            writeCommented(file, "Time/diameter");
            break;
        }
    }

    switch (functionType_)
    {
        case ftMom:
        {
            for (label i = 0; i <= momentOrder_; i++)
            {
                file() << tab << i;
            }

            break;
        }
        default:
        {
            forAll(popBal_.sizeGroups(), sizeGroupi)
            {
                const diameterModels::sizeGroup& fi =
                    *popBal_.sizeGroups()[sizeGroupi];

                switch (abszissaType_)
                {
                    case atDiameter:
                    {
                        file() << tab  << fi.d().value();

                        break;
                    }

                    case atVolume:
                    {
                        file() << tab  << fi.x().value();

                        break;
                    }
                }
            }

            break;
        }
    }

    file << endl;
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
    logFiles(obr_, name),
    dict_(dict),
    regionType_(regionTypeNames_.read(dict.lookup("regionType"))),
    regionName_(word::null),
    functionType_(functionTypeNames_.read(dict.lookup("functionType"))),
    abszissaType_(abszissaTypeNames_.read(dict.lookup("abszissaType"))),
    nCells_(0),
    cellId_(),
    volume_(0.0),
    writeVolume_(dict.lookupOrDefault("writeVolume", false)),
    popBal_
    (
        obr_.lookupObject<Foam::diameterModels::populationBalanceModel>
        (
            dict.lookup("populationBalance")
        )
    ),
    N_(popBal_.sizeGroups().size()),
    momentOrder_(dict.lookupOrDefault<label>("momentOrder", 0)),
    normalize_(dict.lookupOrDefault("normalize", false)),
    sumN_(0.0),
    sumV_(0.0)
{
    read(dict);
    resetName(name);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::sizeDistribution::~sizeDistribution()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::sizeDistribution::read(const dictionary& dict)
{
    if (dict != dict_)
    {
        dict_ = dict;
    }

    fvMeshFunctionObject::read(dict);

    initialise(dict);

    return true;
}


bool Foam::functionObjects::sizeDistribution::execute()
{
    return true;
}


bool Foam::functionObjects::sizeDistribution::write()
{
    logFiles::write();

    if (Pstream::master())
    {
        writeTime(file());
    }

    Log << type() << " " << name() << " write" << nl;

    scalarField V(filterField(mesh().V()));
    combineFields(V);

    sumN_ *= 0.0;
    sumV_ *= 0.0;

    forAll(N_, i)
    {
        const Foam::diameterModels::sizeGroup& fi =
            *popBal_.sizeGroups()[i];

        const volScalarField& alpha = fi.VelocityGroup().phase();

        scalarField Ni(fi*alpha/fi.x());
        scalarField values(filterField(Ni));
        scalarField V(filterField(mesh().V()));

        // Combine onto master
        combineFields(values);
        combineFields(V);

        if (Pstream::master())
        {
            // Calculate volume-averaged number concentration
            N_[i] = sum(V*values)/sum(V);
        }

        sumN_ += N_[i];

        sumV_ += N_[i]*fi.x().value();
    }

    if (Pstream::master())
    {
        switch (functionType_)
        {
            case ftMom:
            {
                for (label m = 0; m <= momentOrder_; m++)
                {
                    scalar result(0.0);

                    forAll(N_, i)
                    {
                        const Foam::diameterModels::sizeGroup& fi =
                            *popBal_.sizeGroups()[i];

                        switch (abszissaType_)
                        {
                            case atVolume:
                            {
                                result += pow(fi.x().value(), m)*N_[i];

                                break;
                            }

                            case atDiameter:
                            {
                                result += pow(fi.d().value(), m)*N_[i];

                                break;
                            }
                        }
                    }

                    file() << tab << result;
                }

                break;
            }

            default:
            {
                forAll(popBal_.sizeGroups(), i)
                {
                    const Foam::diameterModels::sizeGroup& fi =
                        *popBal_.sizeGroups()[i];

                    scalar result(0.0);
                    scalar delta(0.0);

                    switch (abszissaType_)
                    {
                        case atVolume:
                        {
                            delta = popBal_.v()[i+1].value()
                              - popBal_.v()[i].value();

                            break;
                        }

                        case atDiameter:
                        {
                            const scalar& formFactor =
                                fi.VelocityGroup().formFactor().value();

                            delta =
                                pow
                                (
                                    popBal_.v()[i+1].value()
                                   /formFactor,
                                    1.0/3.0
                                )
                              - pow
                                (
                                    popBal_.v()[i].value()
                                   /formFactor,
                                    1.0/3.0
                                );

                            break;
                        }
                    }

                    switch (functionType_)
                    {
                        case ftNdf:
                        {
                            if (normalize_ == true)
                            {
                                result = N_[i]/delta/sumN_;
                            }
                            else
                            {
                                result = N_[i]/delta;
                            }

                            break;
                        }

                        case ftVdf:
                        {
                            if (normalize_ == true)
                            {
                                result = N_[i]*fi.x().value()/delta/sumV_;
                            }
                            else
                            {
                                result = N_[i]*fi.x().value()/delta;
                            }

                            break;
                        }

                        case ftNc:
                        {
                            if (normalize_ == true)
                            {
                                result = N_[i]/sumN_;
                            }
                            else
                            {
                                result = N_[i];
                            }

                            break;
                        }

                        default:
                        {
                            break;
                        }
                    }

                    file()<< tab << result;
                }
            }
        }
    }

    if (Pstream::master())
    {
        file()<< endl;
    }

    Log << endl;

    return true;
}


// ************************************************************************* //
