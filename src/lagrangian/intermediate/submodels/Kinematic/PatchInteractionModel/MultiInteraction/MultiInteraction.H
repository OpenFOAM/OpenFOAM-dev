/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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
    Foam::MultiInteraction

Description
    Runs multiple patch interaction models in turn. Takes dictionary
    where all the subdictionaries are the interaction models.

        // Exit upon first successful interaction or continue doing other
        // models. Returned nteraction status will be true if there has been any
        // interaction (so logical or)
        oneInteractionOnly true;

        model1
        {
            patchInteractionModel coincidentBaffleInteraction;
            coincidentBaffleInteractionCoeffs
            {
                coincidentPatches
                (
                    (pipetteWall_A pipetteCyclic_half0)
                    (pipetteWall_B pipetteCyclic_half1)
                );
            }
        }
        model2
        {
            patchInteractionModel localInteraction;
            localInteractionCoeffs
            {
                patches
                (
                    cWall
                    {
                        type rebound;
                    }
                    pipetteWall_A
                    {
                        type rebound;
                    }
                    pipetteWall_B
                    {
                        type rebound;
                    }
                );
            }
        }


\*---------------------------------------------------------------------------*/

#ifndef MultiInteraction_H
#define MultiInteraction_H

#include "PatchInteractionModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                      Class MultiInteraction Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class MultiInteraction
:
    public PatchInteractionModel<CloudType>
{
    // Private data

        Switch oneInteractionOnly_;

        //- Submodels
        PtrList<PatchInteractionModel<CloudType> > models_;


    // Private Member Functions

        //- Read settings
        bool read(const dictionary&);

public:

    //- Runtime type information
    TypeName("multiInteraction");


    // Constructors

        //- Construct from dictionary
        MultiInteraction(const dictionary& dict, CloudType& cloud);

        //- Construct copy from owner cloud and patch interaction model
        MultiInteraction(const MultiInteraction<CloudType>& pim);

        //- Construct and return a clone using supplied owner cloud
        virtual autoPtr<PatchInteractionModel<CloudType> > clone() const
        {
            return autoPtr<PatchInteractionModel<CloudType> >
            (
                new MultiInteraction<CloudType>(*this)
            );
        }


    //- Destructor
    virtual ~MultiInteraction();


    // Member Functions

        //- Flag to indicate whether model activates patch interaction model
        virtual bool active() const;

        //- Apply velocity correction
        //  Returns true if particle remains in same cell
        virtual bool correct
        (
            typename CloudType::parcelType& p,
            const polyPatch& pp,
            bool& keepParticle,
            const scalar trackFraction,
            const tetIndices& tetIs
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "MultiInteraction.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
