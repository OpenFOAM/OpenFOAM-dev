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
    Foam::functionObjects::patchCutLayerAverage

Description
    This function object writes graphs of patch face values, area-averaged in
    planes perpendicular to a given direction. It adaptively grades the
    distribution of graph points to match the resolution of the mesh.

    Example of function object specification:
    \verbatim
    patchCutLayerAverage1
    {
        type            patchCutLayerAverage;
        libs            ("libpatchCutLayerAverageFunctionObject.so");

        writeControl    writeTime;
        writeInterval   1;

        patch           lowerWall;
        direction       (1 0 0);
        nPoints         100;
        interpolate     no;

        fields          (p U);

        axis            x;
        setFormat       raw;
    }
    \endverbatim

Usage
    \table
        Property    | Description                  | Required | Default value
        type        | type name: patchCutLayerAverage | yes   |
        patch       | Name of the patch            | yes      |
        direction   | Axis along which to plot     | yes      |
        nPoints     | Number of points in the plot | yes      |
        interpolate | Do linear interpolation      | no       | false
        fields      | Fields to plot values of     | yes      |
        axis        | Component of the position to plot against | yes |
        setFormat   | Format of the output file    | yes      |
    \endtable

SourceFiles
    patchCutLayerAverage.C

\*---------------------------------------------------------------------------*/

#ifndef patchCutLayerAverage_H
#define patchCutLayerAverage_H

#include "fvMeshFunctionObject.H"
#include "vector.H"
#include "setWriter.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class polyPatch;

namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                       Class patchCutLayerAverage Declaration
\*---------------------------------------------------------------------------*/

class patchCutLayerAverage
:
    public fvMeshFunctionObject
{
    // Private Structures

        //- Weight for a given face and layer
        struct weight
        {
            label facei, layeri;
            scalar value;
        };


    // Private Data

        //- Name of the patch
        word patchName_;

        //- Direction along which to plot
        vector direction_;

        //- Number of layers
        label nLayers_;

        //- Whether or not to interpolate
        bool interpolate_;

        //- Fields to plot
        wordList fields_;

        //- The direction over which to plot the results
        coordSet::axisType axis_;

        //- File writer
        autoPtr<setWriter> formatter_;

        //- Number of optimisation iterations. Default is 2.
        label nOptimiseIter_;

        //- Weights
        List<weight> weights_;

        //- Layer distances. The distance of the layers' cut-planes to the
        //  origin, in the direction of the given vector. If interpolate is
        //  false, then the distance is taken to the middle of the layer.
        scalarField layerDistances_;

        //- Layer thicknesses. If interpolate is false, then this is the
        //  thickness of the layer that the plot point represents.
        autoPtr<scalarField> layerThicknesses_;

        //- Layer positions. The average position of the layer.
        pointField layerPositions_;


    // Private Member functions

        //- Calculate non-interpolating weights with which to construct plot
        //  values
        List<weight> calcNonInterpolatingWeights
        (
            const scalarField& pointXs,
            const scalarField& faceMinXs,
            const scalarField& faceMaxXs,
            const labelList& faceMinOrder,
            const scalarField& plotXs,
            const bool normalise = true
        ) const;

        //- Calculate interpolating weights with which to construct plot values
        List<weight> calcInterpolatingWeights
        (
            const scalarField& pointXs,
            const scalarField& faceMinXs,
            const scalarField& faceMaxXs,
            const labelList& faceMinOrder,
            const scalarField& plotXs,
            const bool normalise = true
        ) const;

        //- Calculate weights with which to construct plot values
        List<weight> calcWeights
        (
            const scalarField& pointXs,
            const scalarField& faceMinXs,
            const scalarField& faceMaxXs,
            const labelList& faceMinOrder,
            const scalarField& plotXs,
            const bool normalise = true
        ) const;

        //- Construct plot values from face values given a set of weights
        template<class Type>
        inline tmp<Field<Type>> applyWeights
        (
            const List<weight>& weights,
            const Field<Type>& faceValues
        ) const;

        //- Initialise the cached weights
        void initialise();

        //- Return the output path
        fileName outputPath() const;


public:

    //- Runtime type information
    TypeName("patchCutLayerAverage");


    // Constructors

        //- Construct from Time and dictionary
        patchCutLayerAverage
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );


    //- Destructor
    virtual ~patchCutLayerAverage();


    // Member Functions

        //- Read the patchCutLayerAverage data
        virtual bool read(const dictionary&);

        //- Return the list of fields required
        virtual wordList fields() const;

        //- Execute, currently does nothing
        virtual bool execute();

        //- Write the patchCutLayerAverage
        virtual bool write();

        //- Update for mesh point-motion
        virtual void movePoints(const polyMesh&);

        //- Update topology using the given map
        virtual void topoChange(const polyTopoChangeMap&);

        //- Update from another mesh using the given map
        virtual void mapMesh(const polyMeshMap&);

        //- Redistribute or update using the given distribution map
        virtual void distribute(const polyDistributionMap&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
