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
    Foam::VoFphase

Description
    Single incompressible VoFphase derived from the phase-fraction.
    Used as part of the multiPhaseMixture for interface-capturing multi-VoFphase
    simulations.

SourceFiles
    VoFphase.C

\*---------------------------------------------------------------------------*/

#ifndef VoFphase_H
#define VoFphase_H

#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                            Class VoFphase Declaration
\*---------------------------------------------------------------------------*/

class VoFphase
:
    public volScalarField
{
    // Private Data

        word name_;


public:

    // Constructors

        //- Construct from components
        VoFphase
        (
            const word& name,
            const fvMesh& mesh
        );

        //- Return clone
        virtual autoPtr<VoFphase> clone() const;

        //- Return a pointer to a new VoFphase
        //  created on freestore
        //  from Istream
        class iNew
        {
        protected:

            const fvMesh& mesh_;

        public:

            iNew(const fvMesh& mesh)
            :
                mesh_(mesh)
            {}

            virtual ~iNew()
            {}

            virtual autoPtr<VoFphase> operator()(Istream& is) const = 0;
        };


    // Member Functions

        const word& name() const
        {
            return name_;
        }

        const word& keyword() const
        {
            return name();
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
