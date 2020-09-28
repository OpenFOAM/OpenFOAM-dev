/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2013-2020 OpenFOAM Foundation
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

#include "regionSizeDistribution.H"
#include "fvcVolumeIntegrate.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace functionObjects
    {
        defineTypeNameAndDebug(regionSizeDistribution, 0);

        addToRunTimeSelectionTable
        (
            functionObject,
            regionSizeDistribution,
            dictionary
        );
    }

    //- Plus op for FixedList<scalar>
    template<class T, unsigned Size>
    class ListPlusEqOp
    {
        public:
        void operator()
        (
            FixedList<T, Size>& x,
            const FixedList<T, Size>& y
        ) const
        {
            forAll(x, i)
            {
                x[i] += y[i];
            }
        }
    };
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::functionObjects::regionSizeDistribution::writeGraph
(
    const coordSet& coords,
    const word& valueName,
    const scalarField& values
) const
{
    const wordList valNames(1, valueName);

    fileName outputPath = file_.baseTimeDir();
    mkDir(outputPath);

    OFstream str(outputPath/formatterPtr_().getFileName(coords, valNames));

    Info<< "    Writing distribution of " << valueName << " to " << str.name()
        << endl;

    List<const scalarField*> valPtrs(1);
    valPtrs[0] = &values;
    formatterPtr_().write(coords, valNames, valPtrs, str);
}


void Foam::functionObjects::regionSizeDistribution::writeAlphaFields
(
    const regionSplit& regions,
    const Map<label>& patchRegions,
    const Map<scalar>& regionVolume,
    const volScalarField& alpha
) const
{
    const scalar maxDropletVol = 1.0/6.0*pow(maxDiam_, 3);

    // Split alpha field
    // ~~~~~~~~~~~~~~~~~
    // Split into
    //  - liquidCore            : region connected to inlet patches
    //  - per region a volume   : for all other regions
    //  - backgroundAlpha       : remaining alpha


    // Construct field
    volScalarField liquidCore
    (
        IOobject
        (
            alphaName_ + "_liquidCore",
            obr_.time().timeName(),
            obr_,
            IOobject::NO_READ
        ),
        alpha,
        fvPatchField<scalar>::calculatedType()
    );

    volScalarField backgroundAlpha
    (
        IOobject
        (
            alphaName_ + "_background",
            obr_.time().timeName(),
            obr_,
            IOobject::NO_READ
        ),
        alpha,
        fvPatchField<scalar>::calculatedType()
    );


    // Knock out any cell not in patchRegions
    forAll(liquidCore, celli)
    {
        label regionI = regions[celli];
        if (patchRegions.found(regionI))
        {
            backgroundAlpha[celli] = 0;
        }
        else
        {
            liquidCore[celli] = 0;

            scalar regionVol = regionVolume[regionI];
            if (regionVol < maxDropletVol)
            {
                backgroundAlpha[celli] = 0;
            }
        }
    }
    liquidCore.correctBoundaryConditions();
    backgroundAlpha.correctBoundaryConditions();

    Info<< "    Volume of liquid-core = "
        << fvc::domainIntegrate(liquidCore).value()
        << endl;
    Info<< "    Volume of background  = "
        << fvc::domainIntegrate(backgroundAlpha).value()
        << endl;

    Info<< "    Writing liquid-core field to " << liquidCore.name() << endl;
    liquidCore.write();
    Info<< "    Writing background field to " << backgroundAlpha.name() << endl;
    backgroundAlpha.write();
}


Foam::Map<Foam::label>
Foam::functionObjects::regionSizeDistribution::findPatchRegions
(
    const regionSplit& regions
) const
{
    // Mark all regions starting at patches
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    // Count number of patch faces (just for initial sizing)
    const labelHashSet patchIDs(mesh_.boundaryMesh().patchSet(patchNames_));

    label nPatchFaces = 0;
    forAllConstIter(labelHashSet, patchIDs, iter)
    {
        nPatchFaces += mesh_.boundaryMesh()[iter.key()].size();
    }


    Map<label> patchRegions(nPatchFaces);
    forAllConstIter(labelHashSet, patchIDs, iter)
    {
        const polyPatch& pp = mesh_.boundaryMesh()[iter.key()];

        // Collect all regions on the patch
        const labelList& faceCells = pp.faceCells();

        forAll(faceCells, i)
        {
            patchRegions.insert
            (
                regions[faceCells[i]],
                Pstream::myProcNo()     // dummy value
            );
        }
    }


    // Make sure all the processors have the same set of regions
    Pstream::mapCombineGather(patchRegions, minEqOp<label>());
    Pstream::mapCombineScatter(patchRegions);

    return patchRegions;
}


Foam::tmp<Foam::scalarField>
Foam::functionObjects::regionSizeDistribution::divide
(
    const scalarField& num,
    const scalarField& denom
)
{
    tmp<scalarField> tresult(new scalarField(num.size()));
    scalarField& result = tresult.ref();

    forAll(denom, i)
    {
        if (denom[i] != 0)
        {
            result[i] = num[i]/denom[i];
        }
        else
        {
            result[i] = 0.0;
        }
    }
    return tresult;
}


void Foam::functionObjects::regionSizeDistribution::writeGraphs
(
    const word& fieldName,              // name of field
    const labelList& indices,           // index of bin for each region
    const scalarField& sortedField,     // per region field data
    const scalarField& binCount,        // per bin number of regions
    const coordSet& coords              // graph data for bins
) const
{
    if (Pstream::master())
    {
        // Calculate per-bin average
        scalarField binSum(nBins_, 0.0);
        forAll(sortedField, i)
        {
            binSum[indices[i]] += sortedField[i];
        }

        scalarField binAvg(divide(binSum, binCount));

        // Per bin deviation
        scalarField binSqrSum(nBins_, 0.0);
        forAll(sortedField, i)
        {
            binSqrSum[indices[i]] += Foam::sqr(sortedField[i]);
        }
        scalarField binDev
        (
            sqrt(divide(binSqrSum, binCount) - Foam::sqr(binAvg))
        );

        // Write average
        writeGraph(coords, fieldName + "_sum", binSum);
        // Write average
        writeGraph(coords, fieldName + "_avg", binAvg);
        // Write deviation
        writeGraph(coords, fieldName + "_dev", binDev);
    }
}


void Foam::functionObjects::regionSizeDistribution::writeGraphs
(
    const word& fieldName,              // name of field
    const scalarField& cellField,       // per cell field data
    const regionSplit& regions,         // per cell the region(=droplet)
    const labelList& sortedRegions,     // valid regions in sorted order
    const scalarField& sortedNormalisation,

    const labelList& indices,           // per region index of bin
    const scalarField& binCount,        // per bin number of regions
    const coordSet& coords              // graph data for bins
) const
{
    // Sum on a per-region basis. Parallel reduced.
    Map<scalar> regionField(regionSum(regions, cellField));

    // Extract in region order
    scalarField sortedField
    (
        sortedNormalisation
      * extractData
        (
            sortedRegions,
            regionField
        )
    );

    writeGraphs
    (
        fieldName,      // name of field
        indices,        // index of bin for each region
        sortedField,    // per region field data
        binCount,       // per bin number of regions
        coords          // graph data for bins
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::regionSizeDistribution::regionSizeDistribution
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    file_(obr_, name),
    alphaName_(dict.lookup("field")),
    patchNames_(dict.lookup("patches"))
{
    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::regionSizeDistribution::~regionSizeDistribution()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::regionSizeDistribution::read(const dictionary& dict)
{
    dict.lookup("field") >> alphaName_;
    dict.lookup("patches") >> patchNames_;
    dict.lookup("threshold") >> threshold_;
    dict.lookup("maxDiameter") >> maxDiam_;
    minDiam_ = 0.0;
    dict.readIfPresent("minDiameter", minDiam_);
    dict.lookup("nBins") >> nBins_;
    dict.lookup("fields") >> fields_;

    word format(dict.lookup("setFormat"));
    formatterPtr_ = writer<scalar>::New(format);

    if (dict.found("coordinateSystem"))
    {
        coordSysPtr_ = coordinateSystem::New(obr_, dict);

        Info<< "Transforming all vectorFields with coordinate system "
            << coordSysPtr_().name() << endl;
    }

    return true;
}


bool Foam::functionObjects::regionSizeDistribution::execute()
{
    return true;
}


bool Foam::functionObjects::regionSizeDistribution::write()
{
    Info<< type() << " " << name() << " write:" << nl;

    autoPtr<volScalarField> alphaPtr;
    if (obr_.foundObject<volScalarField>(alphaName_))
    {
        Info<< "    Looking up field " << alphaName_ << endl;
    }
    else
    {
        Info<< "    Reading field " << alphaName_ << endl;
        alphaPtr.reset
        (
            new volScalarField
            (
                IOobject
                (
                    alphaName_,
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                ),
                mesh_
            )
        );
    }


    const volScalarField& alpha =
    (
         alphaPtr.valid()
       ? alphaPtr()
       : obr_.lookupObject<volScalarField>(alphaName_)
    );

    Info<< "    Volume of alpha          = "
        << fvc::domainIntegrate(alpha).value()
        << endl;

    const scalar meshVol = gSum(mesh_.V());
    const scalar maxDropletVol = 1.0/6.0*pow(maxDiam_, 3);
    const scalar delta = (maxDiam_-minDiam_)/nBins_;

    Info<< "    Mesh volume              = " << meshVol << endl;
    Info<< "    Maximum droplet diameter = " << maxDiam_ << endl;
    Info<< "    Maximum droplet volume   = " << maxDropletVol << endl;


    // Determine blocked faces
    boolList blockedFace(mesh_.nFaces(), false);
    label nBlocked = 0;

    {
        for (label facei = 0; facei < mesh_.nInternalFaces(); facei++)
        {
            scalar ownVal = alpha[mesh_.faceOwner()[facei]];
            scalar neiVal = alpha[mesh_.faceNeighbour()[facei]];

            if
            (
                (ownVal < threshold_ && neiVal > threshold_)
             || (ownVal > threshold_ && neiVal < threshold_)
            )
            {
                blockedFace[facei] = true;
                nBlocked++;
            }
        }

        // Block coupled faces
        forAll(alpha.boundaryField(), patchi)
        {
            const fvPatchScalarField& fvp = alpha.boundaryField()[patchi];
            if (fvp.coupled())
            {
                tmp<scalarField> townFld(fvp.patchInternalField());
                const scalarField& ownFld = townFld();
                tmp<scalarField> tnbrFld(fvp.patchNeighbourField());
                const scalarField& nbrFld = tnbrFld();

                label start = fvp.patch().patch().start();

                forAll(ownFld, i)
                {
                    scalar ownVal = ownFld[i];
                    scalar neiVal = nbrFld[i];

                    if
                    (
                        (ownVal < threshold_ && neiVal > threshold_)
                     || (ownVal > threshold_ && neiVal < threshold_)
                    )
                    {
                        blockedFace[start+i] = true;
                        nBlocked++;
                    }
                }
            }
        }
    }


    regionSplit regions(mesh_, blockedFace);

    Info<< "    Determined " << regions.nRegions()
        << " disconnected regions" << endl;


    if (debug)
    {
        volScalarField region
        (
            IOobject
            (
                "region",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar(dimless, 0)
        );
        Info<< "    Dumping region as volScalarField to " << region.name()
            << endl;

        forAll(regions, celli)
        {
            region[celli] = regions[celli];
        }
        region.correctBoundaryConditions();
        region.write();
    }


    // Determine regions connected to supplied patches
    Map<label> patchRegions(findPatchRegions(regions));


    // Sum all regions
    const scalarField alphaVol(alpha.primitiveField()*mesh_.V());
    Map<scalar> allRegionVolume(regionSum(regions, mesh_.V()));
    Map<scalar> allRegionAlphaVolume(regionSum(regions, alphaVol));
    Map<label> allRegionNumCells
    (
        regionSum
        (
            regions,
            labelField(mesh_.nCells(), 1.0)
        )
    );

    if (debug)
    {
        Info<< "    " << tab << "Region"
            << tab << "Volume(mesh)"
            << tab << "Volume(" << alpha.name() << "):"
            << tab << "nCells"
            << endl;
        scalar meshSumVol = 0.0;
        scalar alphaSumVol = 0.0;
        label nCells = 0;

        Map<scalar>::const_iterator vIter = allRegionVolume.begin();
        Map<scalar>::const_iterator aIter = allRegionAlphaVolume.begin();
        Map<label>::const_iterator numIter = allRegionNumCells.begin();
        for
        (
            ;
            vIter != allRegionVolume.end()
         && aIter != allRegionAlphaVolume.end();
            ++vIter, ++aIter, ++numIter
        )
        {
            Info<< "    " << tab << vIter.key()
                << tab << vIter()
                << tab << aIter()
                << tab << numIter()
                << endl;

            meshSumVol += vIter();
            alphaSumVol += aIter();
            nCells += numIter();
        }
        Info<< "    " << tab << "Total:"
            << tab << meshSumVol
            << tab << alphaSumVol
            << tab << nCells
            << endl;
        Info<< endl;
    }




    {
        Info<< "    Patch connected regions (liquid core):" << endl;
        Info<< tab << "    Region"
            << tab << "Volume(mesh)"
            << tab << "Volume(" << alpha.name() << "):"
            << endl;
        forAllConstIter(Map<label>, patchRegions, iter)
        {
            label regionI = iter.key();
            Info<< "    " << tab << iter.key()
                << tab << allRegionVolume[regionI]
                << tab << allRegionAlphaVolume[regionI] << endl;

        }
        Info<< endl;
    }

    {
        Info<< "    Background regions:" << endl;
        Info<< "    " << tab << "Region"
            << tab << "Volume(mesh)"
            << tab << "Volume(" << alpha.name() << "):"
            << endl;
        Map<scalar>::const_iterator vIter = allRegionVolume.begin();
        Map<scalar>::const_iterator aIter = allRegionAlphaVolume.begin();

        for
        (
            ;
            vIter != allRegionVolume.end()
         && aIter != allRegionAlphaVolume.end();
            ++vIter, ++aIter
        )
        {
            if
            (
               !patchRegions.found(vIter.key())
             && vIter() >= maxDropletVol
            )
            {
                Info<< "    " << tab << vIter.key()
                    << tab << vIter()
                    << tab << aIter() << endl;
            }
        }
        Info<< endl;
    }



    // Split alpha field
    // ~~~~~~~~~~~~~~~~~
    // Split into
    //  - liquidCore            : region connected to inlet patches
    //  - per region a volume   : for all other regions
    //  - backgroundAlpha       : remaining alpha
    writeAlphaFields(regions, patchRegions, allRegionVolume, alpha);


    // Extract droplet-only allRegionVolume, i.e. delete liquid core
    // (patchRegions) and background regions from maps.
    // Note that we have to use mesh volume (allRegionVolume) and not
    // allRegionAlphaVolume since background might not have alpha in it.
    forAllIter(Map<scalar>, allRegionVolume, vIter)
    {
        label regionI = vIter.key();
        if
        (
            patchRegions.found(regionI)
         || vIter() >= maxDropletVol
        )
        {
            allRegionVolume.erase(vIter);
            allRegionAlphaVolume.erase(regionI);
            allRegionNumCells.erase(regionI);
        }
    }

    if (allRegionVolume.size())
    {
        // Construct mids of bins for plotting
        pointField xBin(nBins_);

        scalar x = 0.5*delta;
        forAll(xBin, i)
        {
            xBin[i] = point(x, 0, 0);
            x += delta;
        }

        const coordSet coords("diameter", "x", xBin, mag(xBin));


        // Get in region order the alpha*volume and diameter
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        const labelList sortedRegions = allRegionAlphaVolume.sortedToc();

        scalarField sortedVols
        (
            extractData
            (
                sortedRegions,
                allRegionAlphaVolume
            )
        );

        // Calculate the diameters
        scalarField sortedDiameters(sortedVols.size());
        forAll(sortedDiameters, i)
        {
            sortedDiameters[i] = Foam::cbrt
            (
                sortedVols[i]
               *6/constant::mathematical::pi
            );
        }

        // Determine the bin index for all the diameters
        labelList indices(sortedDiameters.size());
        forAll(sortedDiameters, i)
        {
            indices[i] = (sortedDiameters[i]-minDiam_)/delta;
        }

        // Calculate the counts per diameter bin
        scalarField binCount(nBins_, 0.0);
        forAll(sortedDiameters, i)
        {
            binCount[indices[i]] += 1.0;
        }

        // Write counts
        if (Pstream::master())
        {
            writeGraph(coords, "count", binCount);
        }

        // Write to log
        {
            Info<< "    Bins:" << endl;
            Info<< "    " << tab << "Bin"
                << tab << "Min diameter"
                << tab << "Count:"
                << endl;

            scalar diam = 0.0;
            forAll(binCount, binI)
            {
                Info<< "    " << tab << binI
                    << tab << diam
                    << tab << binCount[binI] << endl;
                diam += delta;
            }
            Info<< endl;
        }


        // Write average and deviation of droplet volume.
        writeGraphs
        (
            "volume",           // name of field
            indices,            // per region the bin index
            sortedVols,         // per region field data
            binCount,           // per bin number of regions
            coords              // graph data for bins
        );

        // Collect some more field
        {
            wordList scalarNames(obr_.names(volScalarField::typeName));
            labelList selected = findStrings(fields_, scalarNames);

            forAll(selected, i)
            {
                const word& fldName = scalarNames[selected[i]];
                Info<< "    Scalar field " << fldName << endl;

                const scalarField& fld = obr_.lookupObject
                <
                    volScalarField
                >(fldName).primitiveField();

                writeGraphs
                (
                    fldName,            // name of field
                    alphaVol*fld,       // per cell field data

                    regions,            // per cell the region(=droplet)
                    sortedRegions,      // valid regions in sorted order
                    1.0/sortedVols,     // per region normalisation

                    indices,            // index of bin for each region
                    binCount,           // per bin number of regions
                    coords              // graph data for bins
                );
            }
        }
        {
            wordList vectorNames(obr_.names(volVectorField::typeName));
            labelList selected = findStrings(fields_, vectorNames);

            forAll(selected, i)
            {
                const word& fldName = vectorNames[selected[i]];
                Info<< "    Vector field " << fldName << endl;

                vectorField fld = obr_.lookupObject
                <
                    volVectorField
                >(fldName).primitiveField();

                if (coordSysPtr_.valid())
                {
                    Info<< "Transforming vector field " << fldName
                        << " with coordinate system "
                        << coordSysPtr_().name()
                        << endl;

                    fld = coordSysPtr_().localVector(fld);
                }


                // Components

                for (direction cmp = 0; cmp < vector::nComponents; cmp++)
                {
                    writeGraphs
                    (
                        fldName + vector::componentNames[cmp],
                        alphaVol*fld.component(cmp),// per cell field data

                        regions,        // per cell the region(=droplet)
                        sortedRegions,  // valid regions in sorted order
                        1.0/sortedVols, // per region normalisation

                        indices,        // index of bin for each region
                        binCount,       // per bin number of regions
                        coords          // graph data for bins
                    );
                }

                // Magnitude
                writeGraphs
                (
                    fldName + "mag",    // name of field
                    alphaVol*mag(fld),  // per cell field data

                    regions,            // per cell the region(=droplet)
                    sortedRegions,      // valid regions in sorted order
                    1.0/sortedVols,     // per region normalisation

                    indices,            // index of bin for each region
                    binCount,           // per bin number of regions
                    coords              // graph data for bins
                );
            }
        }
    }

    return true;
}


// ************************************************************************* //
