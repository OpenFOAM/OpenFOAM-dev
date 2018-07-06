/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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

#include "layerParameters.H"
#include "polyBoundaryMesh.H"
#include "unitConversion.H"
#include "refinementSurfaces.H"
#include "searchableSurfaces.H"
#include "medialAxisMeshMover.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const Foam::scalar Foam::layerParameters::defaultConcaveAngle = 90;


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::scalar Foam::layerParameters::layerExpansionRatio
(
    const label n,
    const scalar totalOverFirst
) const
{
    if (n <= 1)
    {
        return 1;
    }

    const scalar tol = 1e-8;

    if (mag(n - totalOverFirst) < tol)
    {
        return 1;
    }

    const label maxIters = 100;

    // Calculate the bounds of the solution
    scalar minR;
    scalar maxR;

    if (totalOverFirst < n)
    {
        minR = 0;
        maxR = pow(totalOverFirst/n, 1/(n-1));
    }
    else
    {
        minR = pow(totalOverFirst/n, 1/(n-1));
        maxR = totalOverFirst/(n - 1);
    }

    // Starting guess
    scalar r = 0.5*(minR + maxR);

    for (label i = 0; i < maxIters; ++i)
    {
        const scalar prevr = r;
        const scalar fx = pow(r, n) - totalOverFirst*r - (1 - totalOverFirst);
        const scalar dfx = n*pow(r, n - 1) - totalOverFirst;

        r -= fx/dfx;

        if (mag(r - prevr) < tol)
        {
            break;
        }
    }

    return r;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::layerParameters::layerParameters
(
    const dictionary& dict,
    const polyBoundaryMesh& boundaryMesh
)
:
    dict_(dict),
    numLayers_(boundaryMesh.size(), -1),
    relativeSizes_(dict.lookup("relativeSizes")),
    layerSpec_(ILLEGAL),
    firstLayerThickness_(boundaryMesh.size(), -123),
    finalLayerThickness_(boundaryMesh.size(), -123),
    thickness_(boundaryMesh.size(), -123),
    expansionRatio_(boundaryMesh.size(), -123),
    minThickness_
    (
        boundaryMesh.size(),
        readScalar(dict.lookup("minThickness"))
    ),
    featureAngle_(readScalar(dict.lookup("featureAngle"))),
    concaveAngle_
    (
        dict.lookupOrDefault("concaveAngle", defaultConcaveAngle)
    ),
    nGrow_(readLabel(dict.lookup("nGrow"))),
    maxFaceThicknessRatio_
    (
        readScalar(dict.lookup("maxFaceThicknessRatio"))
    ),
    nBufferCellsNoExtrude_
    (
        readLabel(dict.lookup("nBufferCellsNoExtrude"))
    ),
    nLayerIter_(readLabel(dict.lookup("nLayerIter"))),
    nRelaxedIter_(labelMax),
    additionalReporting_(dict.lookupOrDefault("additionalReporting", false)),
    meshShrinker_
    (
        dict.lookupOrDefault
        (
            "meshShrinker",
            medialAxisMeshMover::typeName
        )
    )
{
    // Detect layer specification mode

    label nSpec = 0;

    bool haveFirst = dict.found("firstLayerThickness");
    if (haveFirst)
    {
        firstLayerThickness_ = scalarField
        (
            boundaryMesh.size(),
            readScalar(dict.lookup("firstLayerThickness"))
        );
        nSpec++;
    }
    bool haveFinal = dict.found("finalLayerThickness");
    if (haveFinal)
    {
        finalLayerThickness_ = scalarField
        (
            boundaryMesh.size(),
            readScalar(dict.lookup("finalLayerThickness"))
        );
        nSpec++;
    }
    bool haveTotal = dict.found("thickness");
    if (haveTotal)
    {
        thickness_ = scalarField
        (
            boundaryMesh.size(),
            readScalar(dict.lookup("thickness"))
        );
        nSpec++;
    }
    bool haveExp = dict.found("expansionRatio");
    if (haveExp)
    {
        expansionRatio_ = scalarField
        (
            boundaryMesh.size(),
            readScalar(dict.lookup("expansionRatio"))
        );
        nSpec++;
    }


    if (haveFirst && haveTotal)
    {
        layerSpec_ = FIRST_AND_TOTAL;
        Info<< "Layer thickness specified as first layer and overall thickness."
            << endl;
    }
    else if (haveFirst && haveExp)
    {
        layerSpec_ = FIRST_AND_EXPANSION;
        Info<< "Layer thickness specified as first layer and expansion ratio."
            << endl;
    }
    else if (haveFinal && haveTotal)
    {
        layerSpec_ = FINAL_AND_TOTAL;
        Info<< "Layer thickness specified as final layer and overall thickness."
            << endl;
    }
    else if (haveFinal && haveExp)
    {
        layerSpec_ = FINAL_AND_EXPANSION;
        Info<< "Layer thickness specified as final layer and expansion ratio."
            << endl;
    }
    else if (haveTotal && haveExp)
    {
        layerSpec_ = TOTAL_AND_EXPANSION;
        Info<< "Layer thickness specified as overall thickness"
            << " and expansion ratio." << endl;
    }


    if (layerSpec_ == ILLEGAL || nSpec != 2)
    {
        FatalIOErrorInFunction
        (
            dict
        )   << "Over- or underspecified layer thickness."
            << " Please specify" << nl
            << "    first layer thickness ('firstLayerThickness')"
            << " and overall thickness ('thickness') or" << nl
            << "    first layer thickness ('firstLayerThickness')"
            << " and expansion ratio ('expansionRatio') or" << nl
            << "    final layer thickness ('finalLayerThickness')"
            << " and expansion ratio ('expansionRatio') or" << nl
            << "    final layer thickness ('finalLayerThickness')"
            << " and overall thickness ('thickness') or" << nl
            << "    overall thickness ('thickness')"
            << " and expansion ratio ('expansionRatio'"
            << exit(FatalIOError);
    }


    dict.readIfPresent("nRelaxedIter", nRelaxedIter_);

    if (nLayerIter_ < 0 || nRelaxedIter_ < 0)
    {
        FatalIOErrorInFunction(dict)
            << "Layer iterations should be >= 0." << endl
            << "nLayerIter:" << nLayerIter_
            << " nRelaxedIter:" << nRelaxedIter_
            << exit(FatalIOError);
    }


    const dictionary& layersDict = dict.subDict("layers");

    forAllConstIter(dictionary, layersDict, iter)
    {
        if (iter().isDict())
        {
            const keyType& key = iter().keyword();
            const labelHashSet patchIDs
            (
                boundaryMesh.patchSet(List<wordRe>(1, wordRe(key)))
            );

            if (patchIDs.size() == 0)
            {
                IOWarningInFunction(layersDict)
                    << "Layer specification for " << key
                    << " does not match any patch." << endl
                    << "Valid patches are " << boundaryMesh.names() << endl;
            }
            else
            {
                const dictionary& layerDict = iter().dict();

                forAllConstIter(labelHashSet, patchIDs, patchiter)
                {
                    const label patchi = patchiter.key();

                    numLayers_[patchi] =
                        readLabel(layerDict.lookup("nSurfaceLayers"));

                    switch (layerSpec_)
                    {
                        case FIRST_AND_TOTAL:
                            layerDict.readIfPresent
                            (
                                "firstLayerThickness",
                                firstLayerThickness_[patchi]
                            );
                            layerDict.readIfPresent
                            (
                                "thickness",
                                thickness_[patchi]
                            );
                        break;

                        case FIRST_AND_EXPANSION:
                            layerDict.readIfPresent
                            (
                                "firstLayerThickness",
                                firstLayerThickness_[patchi]
                            );
                            layerDict.readIfPresent
                            (
                                "expansionRatio",
                                expansionRatio_[patchi]
                            );
                        break;

                        case FINAL_AND_TOTAL:
                            layerDict.readIfPresent
                            (
                                "finalLayerThickness",
                                finalLayerThickness_[patchi]
                            );
                            layerDict.readIfPresent
                            (
                                "thickness",
                                thickness_[patchi]
                            );
                        break;

                        case FINAL_AND_EXPANSION:
                            layerDict.readIfPresent
                            (
                                "finalLayerThickness",
                                finalLayerThickness_[patchi]
                            );
                            layerDict.readIfPresent
                            (
                                "expansionRatio",
                                expansionRatio_[patchi]
                            );
                        break;

                        case TOTAL_AND_EXPANSION:
                            layerDict.readIfPresent
                            (
                                "thickness",
                                thickness_[patchi]
                            );
                            layerDict.readIfPresent
                            (
                                "expansionRatio",
                                expansionRatio_[patchi]
                            );
                        break;

                        default:
                            FatalIOErrorInFunction
                            (
                                dict
                            )   << "problem." << exit(FatalIOError);
                        break;
                    }

                    layerDict.readIfPresent
                    (
                        "minThickness",
                        minThickness_[patchi]
                    );
                }
            }
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::layerParameters::layerThickness
(
    const label nLayers,
    const scalar firstLayerThickess,
    const scalar finalLayerThickess,
    const scalar totalThickness,
    const scalar expansionRatio
) const
{
    switch (layerSpec_)
    {
        case FIRST_AND_TOTAL:
        case FINAL_AND_TOTAL:
        case TOTAL_AND_EXPANSION:
        {
            return totalThickness;
        }
        break;

        case FIRST_AND_EXPANSION:
        {
            if (mag(expansionRatio-1) < small)
            {
                return firstLayerThickess * nLayers;
            }
            else
            {
                return firstLayerThickess
                   *(1 - pow(expansionRatio, nLayers))
                   /(1 - expansionRatio);
            }
        }
        break;

        case FINAL_AND_EXPANSION:
        {
            if (mag(expansionRatio-1) < small)
            {
                return finalLayerThickess * nLayers;
            }
            else
            {
                const scalar invExpansion = 1.0/expansionRatio;

                return finalLayerThickess
                   *(1 - pow(invExpansion, nLayers))
                   /(1 - invExpansion);
            }
        }
        break;

        default:
        {
            FatalErrorInFunction
                << exit(FatalError);
            return -vGreat;
        }
    }
}


Foam::scalar Foam::layerParameters::layerExpansionRatio
(
    const label nLayers,
    const scalar firstLayerThickess,
    const scalar finalLayerThickess,
    const scalar totalThickness,
    const scalar expansionRatio
) const
{
    switch (layerSpec_)
    {
        case FIRST_AND_EXPANSION:
        case FINAL_AND_EXPANSION:
        case TOTAL_AND_EXPANSION:
        {
            return expansionRatio;
        }
        break;

        case FIRST_AND_TOTAL:
        {
            return layerExpansionRatio
            (
                nLayers,
                totalThickness/firstLayerThickess
            );
        }
        break;

        case FINAL_AND_TOTAL:
        {
            return
                1.0
               /layerExpansionRatio
                (
                    nLayers,
                    totalThickness/finalLayerThickess
                );
        }
        break;

        default:
        {
            FatalErrorInFunction
                << "Illegal thickness specification" << exit(FatalError);
            return -vGreat;
        }
    }
}


Foam::scalar Foam::layerParameters::firstLayerThickness
(
    const label nLayers,
    const scalar firstLayerThickess,
    const scalar finalLayerThickess,
    const scalar totalThickness,
    const scalar expansionRatio
) const
{
    switch (layerSpec_)
    {
        case FIRST_AND_EXPANSION:
        case FIRST_AND_TOTAL:
        {
            return firstLayerThickess;
        }

        case FINAL_AND_EXPANSION:
        {
            return finalLayerThickess*pow(1.0/expansionRatio, nLayers-1);
        }
        break;

        case FINAL_AND_TOTAL:
        {
            const scalar r = layerExpansionRatio
            (
                nLayers,
                firstLayerThickess,
                finalLayerThickess,
                totalThickness,
                expansionRatio
            );

            return finalLayerThickess/pow(r, nLayers-1);
        }
        break;

        case TOTAL_AND_EXPANSION:
        {
            const scalar r = finalLayerThicknessRatio
            (
                nLayers,
                expansionRatio
            );

            const scalar finalThickness = r*totalThickness;

            return finalThickness/pow(expansionRatio, nLayers-1);
        }
        break;

        default:
        {
            FatalErrorInFunction
                << "Illegal thickness specification" << exit(FatalError);
            return -vGreat;
        }
    }
}


Foam::scalar Foam::layerParameters::finalLayerThicknessRatio
(
    const label nLayers,
    const scalar expansionRatio
) const
{
    if (nLayers > 0)
    {
        if (mag(expansionRatio-1) < small)
        {
            return 1.0/nLayers;
        }
        else
        {
            return
                pow(expansionRatio, nLayers - 1)
               *(1 - expansionRatio)
               /(1 - pow(expansionRatio, nLayers));
        }
    }
    else
    {
        return 0;
    }
}


// ************************************************************************* //
