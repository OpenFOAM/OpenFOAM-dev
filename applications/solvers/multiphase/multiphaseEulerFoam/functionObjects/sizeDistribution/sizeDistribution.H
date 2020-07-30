/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2017-2020 OpenFOAM Foundation
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
    Foam::functionObjects::sizeDistribution

Description
    This function object calculates and outputs volume-averaged information
    about the size distribution of the dispersed phase, such as the number
    density function or its moments. It is designed to be used exclusively with
    the population balance modeling functionality of the multiphaseEulerFoam
    solver. It can be applied to a specific cellZone or the entire domain. The
    function type determines whether the density function and its moments are
    based on the number of dispersed phase elements in a size group or their
    total volume.

    Example of function object specification:
    \verbatim
    numberDensity
    {
        type           sizeDistribution;
        libs ("libmultiphaseEulerFoamFunctionObjects.so");
        ...
        populationBalance   bubbles;
        regionType          cellZone;
        name                zone0;
        functionType        number;
        coordinateType      volume;
        densityFunction     yes;
    }
    \endverbatim

Usage
    \table
        Property     | Description               | Required | Default value
        type         | type name: sizeDistribution | yes    |
        populationBalance | corresponding populationBalance | yes |
        functionType | number/volume/moments/stdDev | yes   |
        coordinateType | used for density/moment calculation | yes |
        normalize    | normalize concentrations  | no       | no
        densityFunction | compute densityFunction | no      | no
        logBased     | use log of coordinate for density | no | no
        maxOrder     | maxim order of moment output | no    | 3
    \endtable

SourceFiles
    sizeDistribution.C

\*---------------------------------------------------------------------------*/

#ifndef sizeDistribution_H
#define sizeDistribution_H

#include "fvMeshFunctionObject.H"
#include "volRegion.H"
#include "logFiles.H"
#include "populationBalanceModel.H"
#include "writeFile.H"
#include "writer.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                      Class sizeDistribution Declaration
\*---------------------------------------------------------------------------*/

class sizeDistribution
:
    public fvMeshFunctionObject,
    public volRegion,
    public logFiles
{
public:

    // Public Data Types

        //- Function type enumeration
        enum functionType
        {
            ftMoments,
            ftStdDev,
            ftNumber,
            ftVolume
        };

        //- Ordinate type names
        static const NamedEnum<functionType, 4> functionTypeNames_;

        //- Coordinate type enumeration
        enum coordinateType
        {
            ctVolume,
            ctArea,
            ctDiameter,
            ctProjectedAreaDiameter
        };

        //- Coordinate type names
        static const NamedEnum<coordinateType, 4> coordinateTypeNames_;


protected:

    // Protected Data

        //- Reference to fvMesh
        const fvMesh& mesh_;

        //- File containing data for all functionTypes except moments
        writeFile file_;

        //- Output formatter, set to raw
        autoPtr<writer<scalar>> scalarFormatter_;

        //- Reference to populationBalanceModel
        const Foam::diameterModels::populationBalanceModel& popBal_;

        //- Function to evaluate
        functionType functionType_;

        //- Abscissa type
        coordinateType coordinateType_;

        //- List of volume-averaged number concentrations
        scalarField N_;

        //- List of volume-averaged surface areas
        scalarField a_;

        //- List of volume-averaged diameters
        scalarField d_;

        //- Normalize number/volume concentrations
        Switch normalize_;

        //- Determines whether density function is calculated
        Switch densityFunction_;

        //- Geometric standard deviation/density function
        Switch geometric_;

        //- Highest moment order
        label maxOrder_;

        //- Bins (representative volume/area/diameter)
        pointField bins_;

        //- Bin component used according to chosen coordinate type
        label binCmpt_;


    // Protected Member Functions

        //- Filter field according to cellIds
        tmp<scalarField> filterField(const scalarField& field) const;

        //- Correct volume averages
        void correctVolAverages();

        //- Write moments
        void writeMoments();

        //- Write standard deviation
        void writeStdDev();

        //- Write distribution
        void writeDistribution();

        //- Output file header information for functionType moments
        virtual void writeFileHeader(const label i);

        //- Output function for all functionType number/volume
        void writeGraph
        (
            const coordSet& coords,
            const word& functionTypeName,
            const scalarField& values
        );


public:

    //- Runtime type information
    TypeName("sizeDistribution");


    // Constructors

        //- Construct from Time and dictionary
        sizeDistribution
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        sizeDistribution(const sizeDistribution&) = delete;


    //- Destructor
    virtual ~sizeDistribution();


    // Member Functions

        //- Read the sizeDistribution data
        virtual bool read(const dictionary&);

        //- Execute, currently does nothing
        virtual bool execute();

        //- Execute at the final time-loop, currently does nothing
        virtual bool end();

        //- Calculate and write the size distribution
        virtual bool write();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const sizeDistribution&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
