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

Class
    Foam::functionObjects::LagrangianDistribution

Description
    Function to generate a plot of the distribution of the values in a
    Lagrangian field

Usage
    \table
        Property      | Description                  | Required? | Default
        Lagrangian    | Name of the Lagrangian mesh  | yes       |
        field         | Field to operate on          | if fields not specified |
        fields        | List of fields to operate on | if field not specified |
        weightField   | Field with which to weight the distribution | no | none
        weightFields  | List of fields with which to \
                        weight the distribution      | no        | none
        nBins         | The number of bins used in the plot | yes |
        setFormat     | Format of the plot file      | yes       |
    \endtable

    Example specification to generate the diameter distribution:
    \verbatim
    LagrangianDistribution1
    {
        type            LagrangianDistribution;
        libs            ("libLagrangianFunctionObjects.so");
        Lagrangian      cloud;
        field           d;
        weightField     number;
        nBins           20;
        setFormat       raw;
    }
    \endverbatim

SourceFiles
    LagrangianDistribution.C

\*---------------------------------------------------------------------------*/

#ifndef LagrangianDistribution_H
#define LagrangianDistribution_H

#include "LagrangianMeshFunctionObject.H"
#include "setWriter.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                       Class LagrangianDistribution Declaration
\*---------------------------------------------------------------------------*/

class LagrangianDistribution
:
    public LagrangianMeshFunctionObject
{
private:

    // Private Data

        //- List of fields
        wordList fields_;

        //- List of weight fields
        wordList weightFields_;

        //- Number of bins
        label nBins_;

        //- File writer
        autoPtr<setWriter> formatter_;


    // Private Member Functions

        //- Non-virtual read
        void readCoeffs(const dictionary& dict);

        //- Multiply the argument by the given weight field. Return whether or
        //  not the weight field was found.
        template<template<class> class GeoField>
        bool multiplyWeight
        (
            const word& weightFieldName,
            scalarField& weight
        ) const;

        //- Write the distribution of the given scalar field
        void writeDistribution
        (
            const scalarField& weight,
            const word& fieldName,
            const scalarField& field
        );

        //- Write the distribution of a field given its name. Return whether or
        //  not the field was found. In the case of higher rank types, multiple
        //  distributions will be calculated and written independently.
        template<template<class> class GeoField, class Type>
        bool writeDistribution
        (
            const scalarField& weight,
            const word& fieldName
        );


public:

    //- Runtime type information
    TypeName("LagrangianDistribution");


    // Constructors

        //- Construct from Time and dictionary
        LagrangianDistribution
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        LagrangianDistribution(const LagrangianDistribution&) = delete;


    //- Destructor
    virtual ~LagrangianDistribution();


    // Member Functions

        //- Read parameters
        virtual bool read(const dictionary&);

        //- Return the list of fields required
        virtual wordList fields() const;

        //- Execute. Does nothing.
        virtual bool execute();

        //- Write the sum
        virtual bool write();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const LagrangianDistribution&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
