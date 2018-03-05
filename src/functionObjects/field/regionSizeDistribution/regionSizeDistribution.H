/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013-2018 OpenFOAM Foundation
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
    Foam::functionObjects::regionSizeDistribution

Group
    grpFieldFunctionObjects

Description
    Creates a size distribution via interrogating a continuous phase fraction
    field.

    Looks up a phase-fraction (alpha) field and splits the mesh into regions
    based on where the field is below the threshold value.  These
    regions ("droplets") can now be analysed.

    Regions:
    - print the regions connected to a user-defined set of patches.
      (in spray calculation these form the liquid core)
    - print the regions with too large volume.  These are the 'background'
      regions.
    - (debug) write regions as a volScalarField
    - (debug) print for all regions the sum of volume and alpha*volume

    Output (volume scalar) fields include:
    - alpha_liquidCore : alpha with outside liquid core set to 0
    - alpha_background : alpha with outside background set to 0.

    %Histogram:
    - determine histogram of diameter (given minDiameter, maxDiameter, nBins)
    - write graph of number of droplets per bin
    - write graph of sum, average and deviation of droplet volume per bin
    - write graph of sum, average and deviation of user-defined fields.  For
      volVectorFields these are those of the 3 components and the magnitude.

    Example of function object specification:
    \verbatim
    regionSizeDistribution1
    {
        type            regionSizeDistribution;
        libs            ("libfieldFunctionObjects.so");
        ...
        field           alpha;
        patches         (inlet);
        threshold       0.4;
        fields          (p U);
        nBins           100;
        maxDiameter     0.5e-4;
        minDiameter     0;
        setFormat       gnuplot;
        coordinateSystem
        {
            type            cartesian;
            origin          (0 0 0);
            e3              (0 1 1);
            e1              (1 0 0);
        }
    }
    \endverbatim

Usage
    \table
        Property     | Description             | Required    | Default value
        type         | type name: regionSizeDistribution |yes|
        field        | phase field to interrogate | yes      |
        patches      | patches from which the liquid core is identified | yes|
        threshold    | phase fraction applied to delimit regions | yes |
        fields       | fields to sample        | yes         |
        nBins        | number of bins for histogram | yes    |
        maxDiameter  | maximum region equivalent diameter | yes |
        minDiameter  | minimum region equivalent diameter | no  | 0
        setFormat    | writing format          | yes         |
        coordinateSystem | transformation for vector fields | no         |
    \endtable

See also
    Foam::functionObject
    Foam::functionObjects::fvMeshFunctionObject
    Foam::functionObjects::writeFile

SourceFiles
    regionSizeDistribution.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_regionSizeDistribution_H
#define functionObjects_regionSizeDistribution_H

#include "fvMeshFunctionObject.H"
#include "writeFile.H"
#include "writer.H"
#include "Map.H"
#include "volFieldsFwd.H"
#include "wordReList.H"
#include "coordinateSystem.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class regionSplit;

namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                   Class regionSizeDistribution Declaration
\*---------------------------------------------------------------------------*/

class regionSizeDistribution
:
    public fvMeshFunctionObject
{
    // Private data

        writeFile file_;

        //- Name of field
        word alphaName_;

        //- Patches to walk from
        wordReList patchNames_;

        //- Clip value
        scalar threshold_;

        //- Maximum droplet diameter
        scalar maxDiam_;

        //- Minimum droplet diameter
        scalar minDiam_;

        //- Number of bins
        label nBins_;

        //- Names of fields to sample on regions
        wordReList fields_;

        //- Output formatter to write
        autoPtr<writer<scalar>> formatterPtr_;

        //- Optional coordinate system
        autoPtr<coordinateSystem> coordSysPtr_;


    // Private Member Functions

        template<class Type>
        Map<Type> regionSum(const regionSplit&, const Field<Type>&) const;

        //- Get data in order
        template<class Type>
        List<Type> extractData(const UList<label>& keys, const Map<Type>&)
        const;

        void writeGraph
        (
            const coordSet& coords,
            const word& valueName,
            const scalarField& values
        ) const;

        //- Write volfields with the parts of alpha which are not
        //  droplets (liquidCore, backGround)
        void writeAlphaFields
        (
            const regionSplit& regions,
            const Map<label>& keepRegions,
            const Map<scalar>& regionVolume,
            const volScalarField& alpha
        ) const;

        //- Mark all regions starting at patches
        Map<label> findPatchRegions(const regionSplit&) const;

        //- Helper: divide if denom != 0
        static tmp<scalarField> divide(const scalarField&, const scalarField&);

        //- Given per-region data calculate per-bin average/deviation and graph
        void writeGraphs
        (
            const word& fieldName,              // name of field
            const labelList& indices,           // index of bin for each region
            const scalarField& sortedField,     // per region field data
            const scalarField& binCount,        // per bin number of regions
            const coordSet& coords              // graph data for bins
        ) const;

        //- Given per-cell data calculate per-bin average/deviation and graph
        void writeGraphs
        (
            const word& fieldName,              // name of field
            const scalarField& cellField,       // per cell field data

            const regionSplit& regions,         // per cell the region(=droplet)
            const labelList& sortedRegions,     // valid regions in sorted order
            const scalarField& sortedNormalisation,

            const labelList& indices,           // index of bin for each region
            const scalarField& binCount,        // per bin number of regions
            const coordSet& coords              // graph data for bins
        ) const;

        //- Disallow default bitwise copy construct
        regionSizeDistribution(const regionSizeDistribution&);

        //- Disallow default bitwise assignment
        void operator=(const regionSizeDistribution&);


public:

    //- Runtime type information
    TypeName("regionSizeDistribution");


    // Constructors

        //- Construct for given objectRegistry and dictionary.
        //  Allow the possibility to load fields from files
        regionSizeDistribution
        (
            const word& name,
            const Time& runTime,
            const dictionary&
        );


    // Destructor

        virtual ~regionSizeDistribution();


    // Member Functions

        //- Read the regionSizeDistribution data
        virtual bool read(const dictionary&);

        //- Do nothing
        virtual bool execute();

        //- Calculate the regionSizeDistribution and write
        virtual bool write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "regionSizeDistributionTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
