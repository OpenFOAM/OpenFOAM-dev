/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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
    Foam::ParticleErosion

Description
    Function object to create a field of eroded volume, Q, on a specified list
    of patches. The volume is calculated by the model of Finnie et al. The
    implementation follows the description given by the review of Yadav et al.

    References:
    \verbatim
        I Finnie, A Levy, D H McFadden
        "Fundamental Mechanisms of the Erosive Wear of Ductile Metals by Solid
        Particles"
        ASTM STP664, 1979, pages 36-58
    \endverbatim

    \verbatim
        I Finnie and D H McFadden
        "On the Velocity Dependence of the Erosion of Ductile Metal by Solid
        Particle at Low Angles of Incidence"
        Wear, 1978, volume 48, pages 181-190
    \endverbatim

    \verbatim
        G Yadav, S Tiwari, A Rajput, R Jatola, M L Jain
        "A Review: Erosion Wear Models"
        Emerging Trends in Mechanical Engineering, 2016, volume 1, pages 150-154
    \endverbatim

Usage
    \table
        Property    | Description                          | Req'd? | Default
        patches     | The patches on which to calculate Q  | yes    |
        p           | Plastic flow stress                  | yes    |
        psi         | Ratio between depth of contact and length of cut | no | 2
        K           | Ratio of normal and tangential force | no | 2
    \endtable

    Example:
    \verbatim
    <functionName>
    {
        type    particleErosion;
        patches (wall1 wall2);
        p       2.7e9;
    }
    \endverbatim

SourceFiles
    ParticleErosion.C

\*---------------------------------------------------------------------------*/

#ifndef ParticleErosion_H
#define ParticleErosion_H

#include "CloudFunctionObject.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class ParticleErosion Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class ParticleErosion
:
    public CloudFunctionObject<CloudType>
{
    // Private Data

        // Typedefs

            //- Convenience typedef for parcel type
            typedef typename CloudType::parcelType parcelType;


        //- Particle eroded volume field
        autoPtr<volScalarField> QPtr_;

        //- List of patch indices to post-process
        labelList patchIDs_;

        //- Plastic flow stress
        scalar p_;

        //- Ratio between depth of contact and length of cut. Default 2.
        scalar psi_;

        //- Ratio of normal and tangential forces. Default 2.
        scalar K_;


protected:

    // Protected Member Functions

        //- Returns local patchi if patch is in patchIds_ list
        label applyToPatch(const label globalPatchi) const;

        //- Write post-processing info
        virtual void write();


public:

    //- Runtime type information
    TypeName("particleErosion");


    // Constructors

        //- Construct from dictionary
        ParticleErosion
        (
            const dictionary& dict,
            CloudType& owner,
            const word& modelName
        );

        //- Construct copy
        ParticleErosion(const ParticleErosion<CloudType>& pe);

        //- Construct and return a clone
        virtual autoPtr<CloudFunctionObject<CloudType>> clone() const
        {
            return autoPtr<CloudFunctionObject<CloudType>>
            (
                new ParticleErosion<CloudType>(*this)
            );
        }


    //- Destructor
    virtual ~ParticleErosion();


    // Member Functions

        // Evaluation

            //- Pre-evolve hook
            virtual void preEvolve();

            //- Post-patch hook
            virtual void postPatch
            (
                const parcelType& p,
                const polyPatch& pp,
                bool& keepParticle
            );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "ParticleErosion.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
