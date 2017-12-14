/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2017 OpenFOAM Foundation
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
    Foam::singleComponentMixture

Description
    Single component mixture

SourceFiles
    singleComponentMixture.C

\*---------------------------------------------------------------------------*/

#ifndef singleComponentMixture_H
#define singleComponentMixture_H

#include "basicSpecieMixture.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class singleComponentMixture Declaration
\*---------------------------------------------------------------------------*/

template<class ThermoType>
class singleComponentMixture
:
    public basicSpecieMixture
{
    // Private data

        //- Thermo model
        ThermoType thermo_;


public:

    //- The type of thermodynamics this mixture is instantiated for
    typedef ThermoType thermoType;


    // Constructors

        //- Construct from dictionary, mesh and phase name
        singleComponentMixture(const dictionary&, const fvMesh&, const word&);


    //- Destructor
    virtual ~singleComponentMixture();


    // Member Functions

        //- Return the instantiated type name
        //  Note: This defines the same name as pureMixture so that
        //  non-reacting and reacting solvers are compatable
        static word typeName()
        {
            return "pureMixture<" + ThermoType::typeName() + '>';
        }

        //- Get the mixture for the given cell
        const ThermoType& cellMixture(const label celli) const
        {
            return thermo_;
        }

        //- Get the mixture for the given patch face
        const ThermoType& patchFaceMixture
        (
            const label patchi,
            const label facei
        ) const
        {
            return thermo_;
        }

        //- Get the volumetric mixture for the given cell
        const ThermoType& cellVolMixture
        (
            const scalar p,
            const scalar T,
            const label celli
        )
        {
            return thermo_;
        }

        //- Get the volumetric mixture for the given patch face
        const ThermoType& patchFaceVolMixture
        (
            const scalar p,
            const scalar T,
            const label patchi,
            const label facei
        ) const
        {
            return thermo_;
        }

        //- Read dictionary
        void read(const dictionary&);

        //- Return thermo based on index
        inline const ThermoType& getLocalThermo(const label speciei) const
        {
            #ifdef FULLDEBUG
            if (speciei != 0)
            {
                FatalErrorInFunction
                    << "Specie index must be zero for a single component "
                    << "mixture" << exit(FatalError);
            }
            #endif
            return thermo_;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "singleComponentMixture.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
