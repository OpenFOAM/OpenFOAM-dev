/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2017-2022 OpenFOAM Foundation
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

Class
    Foam::functionObjects::populationBalanceSizeDistribution

Description
    Writes out the size distribution determined by a population balance model,
    either for the entire domain or a cell zone. Requires solver post-
    processing.

    The following function object specification for example returns the volume-
    based number density function:

    Example of function object specification:
    \verbatim
    numberDensity
    {
        type                populationBalanceSizeDistribution;
        libs                ("libmultiphaseEulerFoamFunctionObjects.so");
        writeControl        writeTime;
        populationBalance   bubbles;
        functionType        numberDensity;
        coordinateType      volume;
        setFormat           raw;
    }
    \endverbatim

Usage
    \table
        Property          | Description                    | Required  | Default
        populationBalance | population balance name        | yes       |
        functionType      | function type                  | yes       |
        coordinateType    | particle property              | yes       |
        allCoordinates    | write all coordinate values    | no        | false
        normalise         | divide by total concentration  | no        | false
        logTransform      | class width based on log of coordinate\\
                                                           | no        | false
        weightType        | weighting in case of field-dependent particle\\
                            properties                     | no\\
                                                           | numberConcentration
        regionType        | cellZone or all                | no        | all
        name              | name of cellZone if required   | no        |
        setFormat         | output format                  | yes       |
    \endtable

See also
    Foam::diameterModels::populationBalanceModel
    Foam::functionObjects::fvMeshFunctionObject
    Foam::functionObjects::volRegion
    Foam::functionObject

SourceFiles
    populationBalanceSizeDistribution.C

\*---------------------------------------------------------------------------*/

#ifndef populationBalanceSizeDistribution_H
#define populationBalanceSizeDistribution_H

#include "fvMeshFunctionObject.H"
#include "volRegion.H"
#include "populationBalanceModel.H"
#include "writeFile.H"
#include "setWriter.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
             Class populationBalanceSizeDistribution Declaration
\*---------------------------------------------------------------------------*/

class populationBalanceSizeDistribution
:
    public fvMeshFunctionObject,
    public volRegion
{
public:

    // Public Data Types

        //- Function type enumeration
        enum functionType
        {
            numberConcentration,
            numberDensity,
            volumeConcentration,
            volumeDensity,
            areaConcentration,
            areaDensity
        };

        //- Function type names
        static const NamedEnum<functionType, 6> functionTypeNames_;

        //- Coordinate type enumeration
        enum coordinateType
        {
            volume,
            area,
            diameter,
            projectedAreaDiameter
        };

        //- Coordinate type names
        static const NamedEnum<coordinateType, 4> coordinateTypeNames_;

        //- Enumeration for the weight types
        enum class weightType
        {
            numberConcentration,
            volumeConcentration,
            areaConcentration,
            cellVolume
        };

        //- Names of the weight types
        static const NamedEnum<weightType, 4> weightTypeNames_;


private:

    // Private Data

        //- Write file
        writeFile file_;

        //- Reference to mesh
        const fvMesh& mesh_;

        //- Reference to population balance
        const Foam::diameterModels::populationBalanceModel& popBal_;

        //- Function type
        functionType functionType_;

        //- Coordinate type
        coordinateType coordinateType_;

        //- Add values for all coordinate types to output
        Switch allCoordinates_;

        //- Normalise result through division by sum
        Switch normalise_;

        //- Log transform
        Switch logTransform_;

        //- Weight types, relevant if particle properties are field dependent
        weightType weightType_;

        //- Set formatter
        autoPtr<setWriter> formatterPtr_;


    // Private Member Functions

        //- Function type symbolic name for shorter file header
        word functionTypeSymbolicName();

        //- Coordinate type symbolic name for shorter file header
        word coordinateTypeSymbolicName(const coordinateType& cType);

        //- Filter a field according to cellIds
        tmp<scalarField> filterField(const scalarField& field) const;

        //- Field averaged coordinate value
        scalar averageCoordinateValue
        (
            const diameterModels::sizeGroup& fi,
            const coordinateType& cType
        );

        //- Weighted average
        scalar weightedAverage
        (
            const scalarField& fld,
            const diameterModels::sizeGroup& fi
        );


public:

    //- Runtime type information
    TypeName("populationBalanceSizeDistribution");


    // Constructors

        //- Construct from Time and dictionary
        populationBalanceSizeDistribution
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        populationBalanceSizeDistribution
        (
            const populationBalanceSizeDistribution&
        ) = delete;


    //- Destructor
    virtual ~populationBalanceSizeDistribution();


    // Member Functions

        //- Return the list of fields required
        virtual wordList fields() const
        {
            return wordList::null();
        }

        //- Read the populationBalanceSizeDistribution data
        virtual bool read(const dictionary&);

        //- Execute, currently does nothing
        virtual bool execute();

        //- Calculate and write the size distribution
        virtual bool write();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const populationBalanceSizeDistribution&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
