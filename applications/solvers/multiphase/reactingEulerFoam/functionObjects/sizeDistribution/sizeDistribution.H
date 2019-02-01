/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2017-2019 OpenFOAM Foundation
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
    This function object calculates and outputs information about the size
    distribution of the dispersed phase, such as the number density function or
    its moments. It is designed to be used exclusively with the population
    balance modeling functionality of the reactingEulerFoam solvers. It can be
    applied to a specific cellZone or the entire domain.

    Example of function object specification:
    \verbatim
    box.all.numberDensity.volume.bubbles
    {
        type sizeDistribution;
        libs ("libreactingEulerFoamFunctionObjects.so");
        writeControl    outputTime;
        writeInterval   1;
        log             true;
        ...
        functionType        numberDensity;
        abszissaType        volume;
        selectionMode       all;
        populationBalanceModel   bubbles;
        normalize           true;
    }
    \endverbatim

Usage
    \table
        Property        | Description           | Required    | Default value
        type            | type name: sizeDistribution | yes |
        functionType    | numberDensity, volumeDensity, numberConcentration,
                          moments | yes |
        abszissaType    | volume, diameter | yes |
        momentOrder     | Write moment up to given order | no | 0
        selectionMode | Evaluate for cellZone or entire mesh | yes |
        cellZone      | Required if selectionMode is cellZone |   |
        populationBalanceModel    | Respective populationBalanceModel | yes |
        normalize | Normalization | no |
    \endtable

See also
    Foam::diameterModels::populationBalanceModel
    Foam::functionObject
    Foam::functionObjects::fvMeshFunctionObject
    Foam::functionObjects::logFiles

SourceFiles
    sizeDistribution.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_sizeDistribution_H
#define functionObjects_sizeDistribution_H

#include "fvMeshFunctionObject.H"
#include "logFiles.H"
#include "populationBalanceModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class fvMesh;

namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                      Class sizeDistribution Declaration
\*---------------------------------------------------------------------------*/

class sizeDistribution
:
    public fvMeshFunctionObject,
    public logFiles
{

public:

    // Public data types

        //- Selection mode type enumeration
        enum selectionModeTypes
        {
            rtCellZone,
            rtAll
        };

        //- Selection mode type names
        static const NamedEnum<selectionModeTypes, 2> selectionModeTypeNames_;


        //- Function type enumeration
        enum functionTypes
        {
            ftNdf,
            ftVdf,
            ftNc,
            ftMom
        };

        //- Function type names
        static const NamedEnum<functionTypes, 4> functionTypeNames_;


        //- abszissa type enumeration
        enum abszissaTypes
        {
            atDiameter,
            atVolume,
        };

        //- Abszissa type names
        static const NamedEnum<abszissaTypes, 2> abszissaTypeNames_;


protected:

    // Protected data

        //- Construction dictionary
        dictionary dict_;

        //- Selection mode type
        selectionModeTypes selectionModeType_;

        //- Name of selection
        word selectionModeTypeName_;

        //- Function type
        functionTypes functionType_;

        //- Abszissa type
        abszissaTypes abszissaType_;

        //- Global number of cells
        label nCells_;

        //- Local list of cell IDs
        labelList cellId_;

        //- Total volume of the evaluated selection
        scalar volume_;

        //- Optionally write the volume of the sizeDistribution
        bool writeVolume_;

        //- PopulationBalance
        const Foam::diameterModels::populationBalanceModel& popBal_;

        //- Number concentrations
        List<scalar> N_;

        //- Write moments up to specified order with respect to abszissaType
        label momentOrder_;

        //- Normalization switch
        const Switch normalize_;

        //- Sum of number concentrations
        scalar sumN_;

        //- Volumertic sum
        scalar sumV_;


    // Protected Member Functions

        //- Initialise, e.g. cell addressing
        void initialise(const dictionary& dict);

        //- Set cells to evaluate based on a cell zone
        void setCellZoneCells();

        //- Calculate and return volume of the evaluated cell zone
        scalar volume() const;

        //- Combine fields from all processor domains into single field
        void combineFields(scalarField& field);

        //- Filter field according to cellIds
        tmp<scalarField> filterField(const scalarField& field) const;

        //- Output file header information
        virtual void writeFileHeader(const label i);


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


    //- Destructor
    virtual ~sizeDistribution();


    // Member Functions

        //- Return the reference to the construction dictionary
        const dictionary& dict() const
        {
            return dict_;
        }

        //- Return the local list of cell IDs
        const labelList& cellId() const
        {
            return cellId_;
        }

        //- Helper function to return the reference to the mesh
        const fvMesh& mesh() const
        {
            return refCast<const fvMesh>(obr_);
        }

        //- Read from dictionary
        virtual bool read(const dictionary& dict);

        //- Execute
        virtual bool execute();

        //- Write
        virtual bool write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
