/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2023 OpenFOAM Foundation
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
    Foam::SizeDistribution

Description
    Creates graphs of a cloud's size distribution

    Example usage:
    \verbatim
    massFlux1
    {
        type        sizeDistribution;
        nPoints     40;
        setFormat   raw;
    }
    \endverbatim

SourceFiles
    SizeDistribution.C

\*---------------------------------------------------------------------------*/

#ifndef SizeDistribution_H
#define SizeDistribution_H

#include "CloudFunctionObject.H"
#include "setWriter.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class SizeDistribution Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class SizeDistribution
:
    public CloudFunctionObject<CloudType>
{
    // Private Data

        //- Number of points to plot
        const label nPoints_;

        //- Set formatter
        autoPtr<setWriter> formatter_;


protected:

    // Protected Member Functions

        //- Write post-processing info
        virtual void write();


public:

    //- Runtime type information
    TypeName("sizeDistribution");


    // Constructors

        //- Construct from dictionary
        SizeDistribution
        (
            const dictionary& dict,
            CloudType& owner,
            const word& modelName
        );

        //- Construct copy
        SizeDistribution(const SizeDistribution<CloudType>& vf);

        //- Construct and return a clone
        virtual autoPtr<CloudFunctionObject<CloudType>> clone() const
        {
            return autoPtr<CloudFunctionObject<CloudType>>
            (
                new SizeDistribution<CloudType>(*this)
            );
        }


    //- Destructor
    virtual ~SizeDistribution();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "SizeDistribution.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
