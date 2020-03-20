/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2020 OpenFOAM Foundation
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
    Foam::heatTransferPhaseSystem

SourceFiles
    heatTransferPhaseSystem.C

\*---------------------------------------------------------------------------*/

#ifndef heatTransferPhaseSystem_H
#define heatTransferPhaseSystem_H

#include "NamedEnum.H"
#include "primitiveFieldsFwd.H"
#include "volFieldsFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class phasePair;

/*---------------------------------------------------------------------------*\
                   Class heatTransferPhaseSystem Declaration
\*---------------------------------------------------------------------------*/

class heatTransferPhaseSystem
{
public:

    //- Enumeration for the latent heat scheme
    enum class latentHeatScheme
    {
        symmetric,
        upwind
    };

    //- Names of the latent heat schemes
    static const NamedEnum<latentHeatScheme, 2> latentHeatSchemeNames_;


    //- Runtime type information
    TypeName("heatTransferPhaseSystem");


    // Constructors

        //- Default constructor
        heatTransferPhaseSystem();


    //- Destructor
    virtual ~heatTransferPhaseSystem();


    // Member Functions

        //- Return the latent heat for a given pair, mass transfer rate (used
        //  only for it's sign), and interface temperature
        virtual tmp<volScalarField> L
        (
            const phasePair& pair,
            const volScalarField& dmdtf,
            const volScalarField& Tf,
            const latentHeatScheme scheme
        ) const = 0;

        //- As above, but for a cell-set
        virtual tmp<scalarField> L
        (
            const phasePair& pair,
            const scalarField& dmdtf,
            const scalarField& Tf,
            const labelUList& cells,
            const latentHeatScheme scheme
        ) const = 0;

        //- Return the latent heat for a given pair, specie, mass transfer rate
        //  (used only for it's sign), and interface temperature
        virtual tmp<volScalarField> Li
        (
            const phasePair& pair,
            const word& member,
            const volScalarField& dmidtf,
            const volScalarField& Tf,
            const latentHeatScheme scheme
        ) const = 0;

        //- As above, but for a cell-set
        virtual tmp<scalarField> Li
        (
            const phasePair& pair,
            const word& member,
            const scalarField& dmidtf,
            const scalarField& Tf,
            const labelUList& cells,
            const latentHeatScheme scheme
        ) const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
