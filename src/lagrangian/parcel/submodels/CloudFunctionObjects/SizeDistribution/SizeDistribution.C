/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2023-2025 OpenFOAM Foundation
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

#include "SizeDistribution.H"
#include "OSspecific.H"
#include "setWriter.H"

// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

template<class CloudType>
void Foam::SizeDistribution<CloudType>::write()
{
    // Check that there are some parcels
    const label nParcels =
        returnReduce(this->owner().nParcels(), sumOp<label>());
    if (nParcels == 0)
    {
        Info<< type() << ": Not writing the distribution as the cloud "
            << "is empty" << nl << endl;
        return;
    }

    // Check that there is a non-zero range of diameters
    const scalar d0 = this->owner().Dmin(), d1 = this->owner().Dmax();
    if (d1 == d0)
    {
        Info<< type() << ": Not writing the distribution as the cloud "
            << "has uniform particle diameters" << nl << endl;
        return;
    }

    // The x-axis is linearly spaced between the limiting diameters
    scalarField ds(nPoints_);
    forAll(ds, i)
    {
        const scalar f = scalar(i)/(nPoints_ - 1);
        ds[i] = (1 - f)*d0 + f*d1;
    }

    // Calculate the distribution
    scalarField particlePDF(nPoints_, 0), parcelPDF(nPoints_, 0);
    forAllConstIter(typename CloudType, this->owner(), iter)
    {
        const scalar nParticle = iter().nParticle();
        const scalar d = iter().d();

        const scalar f = (d - d0)/(d1 - d0);
        const label i = min(floor(f*(nPoints_ - 1)), nPoints_ - 2);
        const scalar g = f*(nPoints_ - 1) - scalar(i);

        particlePDF[i] += nParticle*(1 - g);
        particlePDF[i + 1] += nParticle*g;

        parcelPDF[i] += 1 - g;
        parcelPDF[i + 1] += g;
    }

    Pstream::listCombineGather(particlePDF, plusEqOp<scalar>());
    Pstream::listCombineScatter(particlePDF);

    Pstream::listCombineGather(parcelPDF, plusEqOp<scalar>());
    Pstream::listCombineScatter(parcelPDF);

    particlePDF.first() *= 2;
    particlePDF.last() *= 2;
    particlePDF /= sum(particlePDF)*(d1 - d0)/(nPoints_ - 1);

    parcelPDF.first() *= 2;
    parcelPDF.last() *= 2;
    parcelPDF /= sum(parcelPDF)*(d1 - d0)/(nPoints_ - 1);

    Info<< type() << ": Writing the distribution to "
        << this->writeTimeDir().relativePath() << nl << endl;

    // Write
    if (Pstream::master())
    {
        mkDir(this->writeTimeDir());

        formatter_->write
        (
            this->writeTimeDir(),
            "distribution",
            coordSet(true, "d", ds),
            "particle-PDF",
            particlePDF,
            "parcel-PDF",
            parcelPDF
        );
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::SizeDistribution<CloudType>::SizeDistribution
(
    const dictionary& dict,
    CloudType& owner,
    const word& modelName
)
:
    CloudFunctionObject<CloudType>(dict, owner, modelName, typeName),
    nPoints_(dict.lookup<label>("nPoints")),
    formatter_(setWriter::New(dict.lookup("setFormat"), dict))
{}


template<class CloudType>
Foam::SizeDistribution<CloudType>::SizeDistribution
(
    const SizeDistribution<CloudType>& vf
)
:
    CloudFunctionObject<CloudType>(vf),
    nPoints_(vf.nPoints_),
    formatter_(vf.formatter_, false)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::SizeDistribution<CloudType>::~SizeDistribution()
{}


// ************************************************************************* //
