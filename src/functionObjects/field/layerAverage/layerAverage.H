/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2021 OpenFOAM Foundation
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
    Foam::layerAverage

Description
    Generates plots of fields averaged over the layers in the mesh

    Example of function object specification:
    \verbatim
    layerAverage1
    {
        type            layerAverage;
        libs            ("libfieldFunctionObjects.so");

        writeControl    writeTime;

        setFormat       raw;

        patches         (bottom);
        zones           (quarterPlane threeQuartersPlane);

        component       y;

        symmetric       true;

        fields          (pMean pPrime2Mean UMean UPrime2Mean k);
    }
    \endverbatim

Usage
    \table
        Property      | Description                | Required   | Default value
        type          | Type name: layerAverage    | yes        |
        setFormat     | Output plotting format     | yes        |
        patches       | Patches that layers extrude from | no | ()
        zones         | Face zones that the layers extrude from | no | ()
        component     | Component of position to plot against | yes |
        symmetric     | Is the geometry symmetric around the centre layer? \
                                                   | no         | false
        field         | Fields to average and plot | yes        |
    \endtable

SourceFiles
    layerAverage.C
    layerAverageTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef layerAverage_H
#define layerAverage_H

#include "fvMeshFunctionObject.H"
#include "setWriter.H"
#include "boolList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                          Class layerAverage Declaration
\*---------------------------------------------------------------------------*/

class layerAverage
:
    public fvMeshFunctionObject
{
    // Private Static Data

        //- Names of vector components
        static const NamedEnum<vector::components, 3> axisNames_;


    // Private Data

        //- Patches which form the start of the layers
        labelList patchIDs_;

        //- Zones which form the start of the layers
        labelList zoneIDs_;

        //- Zones on which the layers are considered to end
        labelList endZoneIDs_;

        //- Is the case symmetric?
        bool symmetric_;

        //- The direction over which to plot the results
        vector::components axis_;

        //- Per cell the global layer
        label nLayers_;

        //- Per cell the global layer
        labelList cellLayer_;

        //- Per global layer the number of cells
        scalarField layerCount_;

        //- From sorted layer back to unsorted global layer
        labelList sortMap_;

        //- Sorted component of cell centres
        scalarField x_;

        //- Fields to sample
        wordList fields_;

        //- Set formatter
        autoPtr<setWriter> formatter_;


    // Private Member Functions

        //- Walk through layers marking faces that separate layers and cells
        //  that are within layers
        void walkOppositeFaces
        (
            const labelList& startFaces,
            const boolList& startFaceIntoOwners,
            boolList& blockedFace,
            boolList& cellIsLayer
        );

        //- Create the layer information, the sort map, and the scalar axis
        void calcLayers();

        //- Return the coefficient to multiply onto symmetric values
        template<class T>
        T symmetricCoeff() const;

        //- Sum field per layer
        template<class T>
        Field<T> sum(const Field<T>& cellField) const;

        //- Average a field per layer
        template<class T>
        Field<T> average(const Field<T>& cellField) const;


public:

    //- Runtime type information
    TypeName("layerAverage");


    // Constructors

        //- Construct from Time and dictionary
        layerAverage
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        layerAverage(const layerAverage&) = delete;


    //- Destructor
    virtual ~layerAverage();


    // Member Functions

        //- Read the field average data
        virtual bool read(const dictionary&);

        //- Return the list of fields required
        virtual wordList fields() const;

        //- Do nothing
        virtual bool execute();

        //- Calculate and write the graphs
        virtual bool write();

        //- Update for changes of mesh
        virtual void updateMesh(const mapPolyMesh&);

        //- Update for mesh point-motion
        virtual void movePoints(const polyMesh&);


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const layerAverage&) = delete;
};


template<>
vector layerAverage::symmetricCoeff<vector>() const;

template<>
symmTensor layerAverage::symmetricCoeff<symmTensor>() const;

template<>
tensor layerAverage::symmetricCoeff<tensor>() const;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "layerAverageTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
