/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013 OpenFOAM Foundation
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
#include "volFields.H"
#include "regionSplit.H"
#include "fvcVolumeIntegrate.H"
#include "mathematicalConstants.H"
#include "stringListOps.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(regionSizeDistribution, 0);

    //- plus op for FixedList<scalar>
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

void Foam::regionSizeDistribution::writeGraph
(
    const coordSet& coords,
    const word& valueName,
    const scalarField& values
) const
{
    const wordList valNames(1, valueName);

    fileName outputPath = baseTimeDir();
    mkDir(outputPath);

    OFstream str(outputPath/formatterPtr_().getFileName(coords, valNames));

    Info<< "Writing distribution of " << valueName << " to " << str.name()
        << endl;

    List<const scalarField*> valPtrs(1);
    valPtrs[0] = &values;
    formatterPtr_().write(coords, valNames, valPtrs, str);
}


void Foam::regionSizeDistribution::writeAlphaFields
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
    forAll(liquidCore, cellI)
    {
        label regionI = regions[cellI];
        if (patchRegions.found(regionI))
        {
            backgroundAlpha[cellI] = 0;
        }
        else
        {
            liquidCore[cellI] = 0;

            scalar regionVol = regionVolume[regionI];
            if (regionVol < maxDropletVol)
            {
                backgroundAlpha[cellI] = 0;
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

    Info<< "Writing liquid-core field to " << liquidCore.name() << endl;
    liquidCore.write();
    Info<< "Writing background field to " << backgroundAlpha.name() << endl;
    backgroundAlpha.write();
}


Foam::Map<Foam::label> Foam::regionSizeDistribution::findPatchRegions
(
    const polyMesh& mesh,
    const regionSplit& regions
) const
{
    // Mark all regions starting at patches
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    // Count number of patch faces (just for initial sizing)
    const labelHashSet patchIDs(mesh.boundaryMesh().patchSet(patchNames_));

    label nPatchFaces = 0;
    forAllConstIter(labelHashSet, patchIDs, iter)
    {
        nPatchFaces += mesh.boundaryMesh()[iter.key()].size();
    }


    Map<label> patchRegions(nPatchFaces);
    forAllConstIter(labelHashSet, patchIDs, iter)
    {
        const polyPatch& pp = mesh.boundaryMesh()[iter.key()];

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


Foam::tmp<Foam::scalarField> Foam::regionSizeDistribution::divide
(
    const scalarField& num,
    const scalarField& denom
)
{
    tmp<scalarField> tresult(new scalarField(num.size()));
    scalarField& result = tresult();

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


void Foam::regionSizeDistribution::writeGraphs
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


void Foam::regionSizeDistribution::writeGraphs
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

Foam::regionSizeDistribution::regionSizeDistribution
(
    const word& name,
    const objectRegistry& obr,
    const dictionary& dict,
    const bool loadFromFiles
)
:
    functionObjectFile(obr, name, typeName),
    name_(name),
    obr_(obr),
    active_(true),
    alphaName_(dict.lookup("field")),
    patchNames_(dict.lookup("patches"))
{
    // Check if the available mesh is an fvMesh, otherwise deactivate
    if (isA<fvMesh>(obr_))
    {
        read(dict);
    }
    else
    {
        active_ = false;
        WarningIn
        (
            "regionSizeDistribution::regionSizeDistribution"
            "("
                "const word&,  "
                "const objectRegistry&, "
                "const dictionary&, "
                "const bool"
            ")"
        )   << "No fvMesh available, deactivating " << name_ << nl
            << endl;
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::regionSizeDistribution::~regionSizeDistribution()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::regionSizeDistribution::read(const dictionary& dict)
{
    if (active_)
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
            coordSysPtr_.reset(new coordinateSystem(obr_, dict));

            Info<< "Transforming all vectorFields with coordinate system "
                << coordSysPtr_().name() << endl;
        }
    }
}


void Foam::regionSizeDistribution::execute()
{
    // Do nothing - only valid on write
}


void Foam::regionSizeDistribution::end()
{
    // Do nothing - only valid on write
}


void Foam::regionSizeDistribution::timeSet()
{
    // Do nothing - only valid on write
}


void Foam::regionSizeDistribution::write()
{
    if (active_)
    {
        Info<< type() << " " << name_ << " output:" << nl;

        const fvMesh& mesh = refCast<const fvMesh>(obr_);

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
                        mesh.time().timeName(),
                        mesh,
                        IOobject::MUST_READ,
                        IOobject::NO_WRITE
                    ),
                    mesh
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

        const scalar meshVol = gSum(mesh.V());
        const scalar maxDropletVol = 1.0/6.0*pow(maxDiam_, 3);
        const scalar delta = (maxDiam_-minDiam_)/nBins_;

        Info<< "    Mesh volume              = " << meshVol << endl;
        Info<< "    Maximum droplet diameter = " << maxDiam_ << endl;
        Info<< "    Maximum droplet volume   = " << maxDropletVol << endl;


        // Determine blocked faces
        boolList blockedFace(mesh.nFaces(), false);
        label nBlocked = 0;

        {
            for (label faceI = 0; faceI < mesh.nInternalFaces(); faceI++)
            {
                scalar ownVal = alpha[mesh.faceOwner()[faceI]];
                scalar neiVal = alpha[mesh.faceNeighbour()[faceI]];

                if
                (
                    (ownVal < threshold_ && neiVal > threshold_)
                 || (ownVal > threshold_ && neiVal < threshold_)
                )
                {
                    blockedFace[faceI] = true;
                    nBlocked++;
                }
            }

            // Block coupled faces
            forAll(alpha.boundaryField(), patchI)
            {
                const fvPatchScalarField& fvp = alpha.boundaryField()[patchI];
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


        regionSplit regions(mesh, blockedFace);

        Info<< "    Determined " << regions.nRegions()
            << " disconnected regions" << endl;


        if (debug)
        {
            volScalarField region
            (
                IOobject
                (
                    "region",
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh,
                dimensionedScalar("zero", dimless, 0)
            );
            Info<< "    Dumping region as volScalarField to " << region.name()
                << endl;

            forAll(regions, cellI)
            {
                region[cellI] = regions[cellI];
            }
            region.correctBoundaryConditions();
            region.write();
        }


        // Determine regions connected to supplied patches
        Map<label> patchRegions(findPatchRegions(mesh, regions));



        // Sum all regions
        const scalarField alphaVol(alpha.internalField()*mesh.V());
        Map<scalar> allRegionVolume(regionSum(regions, mesh.V()));
        Map<scalar> allRegionAlphaVolume(regionSum(regions, alphaVol));
        Map<label> allRegionNumCells
        (
            regionSum
            (
                regions,
                labelField(mesh.nCells(), 1.0)
            )
        );

        if (debug)
        {
            Info<< "    " << token::TAB << "Region"
                << token::TAB << "Volume(mesh)"
                << token::TAB << "Volume(" << alpha.name() << "):"
                << token::TAB << "nCells"
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
                Info<< "    " << token::TAB << vIter.key()
                    << token::TAB << vIter()
                    << token::TAB << aIter()
                    << token::TAB << numIter()
                    << endl;

                meshSumVol += vIter();
                alphaSumVol += aIter();
                nCells += numIter();
            }
            Info<< "    " << token::TAB << "Total:"
                << token::TAB << meshSumVol
                << token::TAB << alphaSumVol
                << token::TAB << nCells
                << endl;
            Info<< endl;
        }




        {
            Info<< "    Patch connected regions (liquid core):" << endl;
            Info<< token::TAB << "    Region"
                << token::TAB << "Volume(mesh)"
                << token::TAB << "Volume(" << alpha.name() << "):"
                << endl;
            forAllConstIter(Map<label>, patchRegions, iter)
            {
                label regionI = iter.key();
                Info<< "    " << token::TAB << iter.key()
                    << token::TAB << allRegionVolume[regionI]
                    << token::TAB << allRegionAlphaVolume[regionI] << endl;

            }
            Info<< endl;
        }

        {
            Info<< "    Background regions:" << endl;
            Info<< "    " << token::TAB << "Region"
                << token::TAB << "Volume(mesh)"
                << token::TAB << "Volume(" << alpha.name() << "):"
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
                    Info<< "    " << token::TAB << vIter.key()
                        << token::TAB << vIter()
                        << token::TAB << aIter() << endl;
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

            // Write to screen
            {
                Info<< "    Bins:" << endl;
                Info<< "    " << token::TAB << "Bin"
                    << token::TAB << "Min diameter"
                    << token::TAB << "Count:"
                    << endl;

                scalar diam = 0.0;
                forAll(binCount, binI)
                {
                    Info<< "    " << token::TAB << binI
                        << token::TAB << diam
                        << token::TAB << binCount[binI] << endl;
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
                    >(fldName).internalField();

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
                    >(fldName).internalField();

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
    }
}


// ************************************************************************* //
