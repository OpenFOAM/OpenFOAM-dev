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
    Foam::multiphaseVoFMixture

Description
    Multiphase VoF mixture with support for interface properties

    Surface tension and contact-angle is handled for the interface
    between each phase-pair.

SourceFiles
    multiphaseVoFMixture.C

\*---------------------------------------------------------------------------*/

#ifndef multiphaseVoFMixture_H
#define multiphaseVoFMixture_H

#include "VoFMixture.H"
#include "IOdictionary.H"
#include "PtrListDictionary.H"
#include "VoFphase.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class multiphaseVoFMixture Declaration
\*---------------------------------------------------------------------------*/

class multiphaseVoFMixture
:
    public VoFMixture,
    public IOdictionary
{
public:

    class interfacePair
    :
        public Pair<word>
    {
    public:

        class hash
        :
            public Hash<interfacePair>
        {
        public:

            hash()
            {}

            label operator()(const interfacePair& key) const
            {
                return word::hash()(key.first()) + word::hash()(key.second());
            }
        };


        // Constructors

            interfacePair()
            {}

            interfacePair(const word& alpha1Name, const word& alpha2Name)
            :
                Pair<word>(alpha1Name, alpha2Name)
            {}

            interfacePair(const VoFphase& alpha1, const VoFphase& alpha2)
            :
                Pair<word>(alpha1.name(), alpha2.name())
            {}


        // Friend Operators

            friend bool operator==
            (
                const interfacePair& a,
                const interfacePair& b
            )
            {
                return
                (
                    ((a.first() == b.first()) && (a.second() == b.second()))
                 || ((a.first() == b.second()) && (a.second() == b.first()))
                );
            }

            friend bool operator!=
            (
                const interfacePair& a,
                const interfacePair& b
            )
            {
                return (!(a == b));
            }
    };


protected:

    // Protected Data

        //- Dictionary of phases
        PtrListDictionary<VoFphase> phases_;

        const fvMesh& mesh_;

        volScalarField alphas_;

        typedef HashTable<scalar, interfacePair, interfacePair::hash>
            sigmaTable;

        sigmaTable sigmas_;
        dimensionSet dimSigma_;

        //- Stabilisation for normalisation of the interface normal
        const dimensionedScalar deltaN_;


    // Private Member Functions

        void calcAlphas();

        tmp<surfaceVectorField> nHatfv
        (
            const volScalarField& alpha1,
            const volScalarField& alpha2
        ) const;

        tmp<volScalarField> K
        (
            const VoFphase& alpha1,
            const VoFphase& alpha2,
            const volVectorField& U
        ) const;


public:

    TypeName("multiphaseVoFMixture");


    // Constructors

        //- Construct from fvMesh
        multiphaseVoFMixture
        (
            const fvMesh& mesh,
            const VoFphase::iNew&
        );


    //- Destructor
    virtual ~multiphaseVoFMixture()
    {}


    // Member Functions

        //- Return the phases
        PtrListDictionary<VoFphase>& phases()
        {
            return phases_;
        }

        tmp<surfaceScalarField> nHatf
        (
            const volScalarField& alpha1,
            const volScalarField& alpha2
        ) const;

        tmp<surfaceScalarField> surfaceTensionForce
        (
            const volVectorField& U
        ) const;

        //- Indicator of the proximity of the interface
        //  Field values are 1 near and 0 away for the interface.
        tmp<volScalarField> nearInterface() const;

        //- Read base phaseProperties dictionary
        bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
