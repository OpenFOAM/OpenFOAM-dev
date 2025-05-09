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
    Foam::functionObjects::cloudSurfaceDistribution

Description
    Function to generate a plot of the distribution of the values of particles
    that pass through a face-zone, face-set or patch

Usage
    \table
        Property      | Description                  | Required? | Default
        cloud         | Name of the cloud            | yes       |
        select        | The type of face selection; \
                        faceZone, faceSet or patch   | if none or multiple \
                                                       selections specified |
        faceZone      | The face zone for which the \
                        PDF is calculated            | if select is faceZone |
        faceSet       | The face set for which the \
                        PDF is calculated            | if select is faceSet |
        patch         | The patch for which the \
                        PDF is calculated            | is select is patch |
        field         | Field to operate on          | if fields not specified |
        fields        | List of fields to operate on | if field not specified |
        weightField   | Field with which to weight the distribution | no | none
        weightFields  | List of fields with which to \
                        weight the distribution      | no        | none
        nBins         | The number of bins used in the plot(s) | yes |
        setFormat     | Format of the plot file      | yes       |
    \endtable

    Example specification to generate the diameter distribution:
    \verbatim
    cloudSurfaceDistribution1
    {
        type            cloudSurfaceDistribution;
        libs            ("libLagrangianCloudFunctionObjects.so");
        cloud           cloud;
        patch           outlet;
        field           d;
        weightField     number;
        nBins           20;
        setFormat       raw;
    }
    \endverbatim

SourceFiles
    cloudSurfaceDistribution.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_cloudSurfaceDistribution_H
#define functionObjects_cloudSurfaceDistribution_H

#include "cloudFunctionObject.H"
#include "setWriter.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                       Class cloudSurfaceDistribution Declaration
\*---------------------------------------------------------------------------*/

class cloudSurfaceDistribution
:
    public LagrangianMeshFunctionObject,
    public cloudFunctionObject
{
public:

    // Public Data Types

        //- Selection type enumeration
        enum class selectionType
        {
            faceZone,
            faceSet,
            patch
        };

        //- Selection type names
        static const NamedEnum<selectionType, 3> selectionTypeNames;


private:

    // Private Data

        //- List of fields
        wordList fields_;

        //- List of weight fields
        wordList weightFields_;

        //- Selection type
        selectionType selectionType_;

        //- Name of face selection
        word selectionName_;

        //- Stored set of selected faces. Used if the selection type is faceSet.
        labelHashSet selectionSet_;

        //- Number of bins
        label nBins_;

        //- File writer
        autoPtr<setWriter> formatter_;

        //- Sums for each component of each field
        List<List<scalarField>> sums_;

        //- The ranges of each component of each field
        List<List<Pair<scalar>>> ranges_;

        //- Total number of samples taken for each field
        labelList nSamples_;


    // Private Member Functions

        //- Read a list of field names
        static wordList readFields
        (
            const dictionary& dict,
            const word& key,
            const wordList& defaultValue = NullObjectRef<wordList>()
        );

        //- Read the selection type
        selectionType readSelectionType(const dictionary& dict);

        //- Get the component names for a given field
        const char* const* componentNames(const label fieldi) const;

        //- Non-virtual read
        void readCoeffs(const dictionary& dict, const bool props);

        //- Return the properties dictionary IO
        IOobject propsDictIo(const IOobject::readOption r) const;

        //- Return a list of flags to indicate which particles in the sub-mesh
        //  should be included in the PDF calculation
        boolList selected(const LagrangianSubScalarSubField& fraction) const;

        //- Multiply the given weight field to the current weight product
        template<template<class> class GeoField>
        bool multiplyWeight
        (
            const LagrangianSubMesh& subMesh,
            const label weightFieldi,
            scalarField& weight
        ) const;

        //- Add a field to its PDF
        template<template<class> class GeoField, class Type>
        bool addField
        (
            const LagrangianSubMesh& subMesh,
            const boolList& selected,
            const scalarField& weight,
            const label fieldi
        );

        //- Write the distribution of the given scalar field
        void writeDistribution
        (
            const word& fieldName,
            const word& componentName,
            const scalarField& x,
            const scalarField& PDF,
            const scalarField& CDF
        ) const;


public:

    //- Runtime type information
    TypeName("cloudSurfaceDistribution");


    // Constructors

        //- Construct from Time and dictionary
        cloudSurfaceDistribution
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        cloudSurfaceDistribution(const cloudSurfaceDistribution&) = delete;


    //- Destructor
    virtual ~cloudSurfaceDistribution();


    // Member Functions

        //- Read parameters
        virtual bool read(const dictionary&);

        //- Return the list of fields required
        virtual wordList fields() const;

        //- Return false so this function does not execute at the start
        virtual bool executeAtStart() const;

        //- Do nothing. Everything happens in faces crossing hooks.
        virtual bool execute();

        //- Hook before face crossings of a specific sub-mesh
        virtual void preCrossFaces
        (
            const LagrangianSubScalarSubField& fraction
        );

        //- Hook following face crossings of a specific sub-mesh
        virtual void postCrossFaces
        (
            const LagrangianSubScalarSubField& fraction
        );

        //- Write the number flux
        virtual bool write();

        //- Clear the number flux
        virtual bool clear();

        //- Update topology using the given map
        virtual void topoChange(const polyTopoChangeMap& map);

        //- Update from another mesh using the given map
        virtual void mapMesh(const polyMeshMap&);

        //- Redistribute or update using the given distribution map
        virtual void distribute(const polyDistributionMap&);


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const cloudSurfaceDistribution&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
