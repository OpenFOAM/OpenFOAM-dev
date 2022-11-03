/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2022 OpenFOAM Foundation
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
    Foam::functionObjects::populationBalanceMoments

Description
    Calculates and writes out integral (integer moments) or mean properties
    (mean, variance, standard deviation) of a size distribution determined by a
    population balance model. Requires solver post-processing.

    The following function object specification for example returns the first
    moment of the volume-based number density function which is equivalent to
    the phase fraction of the particulate phase:

    \verbatim
    populationBalanceMoments
    {
        type                populationBalanceMoments;
        libs                ("libmultiphaseEulerFoamFunctionObjects.so");
        executeControl      timeStep;
        writeControl        writeTime;
        populationBalance   bubbles;
        momentType          integerMoment;
        coordinateType      volume;
        order               1;
    }
    \endverbatim

Usage
    \table
        Property          | Description              | Required  | Default
        populationBalance | population balance name  | yes       |
        momentType        | desired moment of the distribution\\
                                                     | yes       |
        coordinateType    | particle property        | yes       |
        weightType        | number/volume/area concentration\\
                                                     | no\\
                                                     | numberConcentration
        order             | order of integer moment  | for integer moments |
        meanType          | arithmetic or geometric  | for non-integer moments\\
                                                                 | arithmetic
    \endtable

See also
    Foam::diameterModels::populationBalanceModel
    Foam::functionObjects::fvMeshFunctionObject
    Foam::functionObject

SourceFiles
    populationBalanceMoments.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_populationBalanceMoments_H
#define functionObjects_populationBalanceMoments_H

#include "fvMeshFunctionObject.H"
#include "populationBalanceModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                  Class populationBalanceMoments Declaration
\*---------------------------------------------------------------------------*/

class populationBalanceMoments
:
    public fvMeshFunctionObject
{
public:

    //- Enumeration for the moment types
    enum class momentType
    {
        integerMoment,
        mean,
        variance,
        stdDev
    };

    //- Names of the moment types
    static const NamedEnum<momentType, 4> momentTypeNames_;

    //- Enumeration for the coordinate types
    enum class coordinateType
    {
        volume,
        area,
        diameter
    };

    //- Names of the coordinate types
    static const NamedEnum<coordinateType, 3> coordinateTypeNames_;

    //- Enumeration for the weight types
    enum class weightType
    {
        numberConcentration,
        volumeConcentration,
        areaConcentration
    };

    //- Names of the weight types
    static const NamedEnum<weightType, 3> weightTypeNames_;

    //- Enumeration for the mean types
    enum class meanType
    {
        arithmetic,
        geometric,
        notApplicable
    };

    //- Names of the mean types
    static const NamedEnum<meanType, 3> meanTypeNames_;


private:

    // Private Data

        //- Reference to population balance
        const Foam::diameterModels::populationBalanceModel& popBal_;

        //- Moment type
        momentType momentType_;

        //- Coordinate type
        coordinateType coordinateType_;

        //- Weight type
        weightType weightType_;

        //- Mean type
        meanType meanType_;

        //- Integer moment order
        int order_;

        //- Result field
        autoPtr<volScalarField> fldPtr_;


    // Private Member Functions

        //- Coordinate type symbolic name for shorter field names
        word coordinateTypeSymbolicName();

        //- Weight type symbolic name for shorter field names
        word weightTypeSymbolicName();

        //- Default field name
        word defaultFldName();

        //- Integer moment field name
        word integerMomentFldName();

        //- Set dimensions
        void setDimensions(volScalarField& fld, momentType momType);

        //- Total concentration
        tmp<volScalarField> totalConcentration();

        //- Mean value
        tmp<volScalarField> mean();

        //- Variance
        tmp<volScalarField> variance();

        //- Standard deviation
        tmp<volScalarField> stdDev();


public:

    //- Runtime type information
    TypeName("populationBalanceMoments");


    // Constructors

        //- Construct from Time and dictionary
        populationBalanceMoments
        (
            const word& name,
            const Time& runTime,
            const dictionary&
        );

        //- Disallow default bitwise copy construction
        populationBalanceMoments(const populationBalanceMoments&) = delete;


    //- Destructor
    virtual ~populationBalanceMoments();


    // Member Functions

        //- Read the data
        virtual bool read(const dictionary&);

        //- Return the list of fields required
        virtual wordList fields() const
        {
            return wordList::null();
        }

        //- Calculate the moment fields
        virtual bool execute();

        //- Write the moment fields
        virtual bool write();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const populationBalanceMoments&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
