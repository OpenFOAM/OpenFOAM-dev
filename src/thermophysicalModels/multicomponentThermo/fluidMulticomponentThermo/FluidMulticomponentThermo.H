/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2023 OpenFOAM Foundation
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
    Foam::FluidMulticomponentThermo

Description
    Fluid multi-component thermo implementation

SourceFiles
    FluidMulticomponentThermo.C

\*---------------------------------------------------------------------------*/

#ifndef FluidMulticomponentThermo_H
#define FluidMulticomponentThermo_H

#include "BasicThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class FluidMulticomponentThermo Declaration
\*---------------------------------------------------------------------------*/

template<class BaseThermo>
class FluidMulticomponentThermo
:
    public BaseThermo
{
public:

    // Constructors

        //- Construct from mesh and phase name
        FluidMulticomponentThermo(const fvMesh&, const word& phaseName);

        //- Disallow default bitwise copy construction
        FluidMulticomponentThermo
        (
            const FluidMulticomponentThermo<BaseThermo>&
        ) = delete;


    //- Destructor
    virtual ~FluidMulticomponentThermo();


    // Member Functions

        // Specie transport properties

            //- Dynamic viscosity [kg/m/s]
            virtual scalar mui
            (
                const label speciei,
                const scalar p,
                const scalar T
            ) const;

            //- Dynamic viscosity [kg/m/s]
            virtual tmp<volScalarField> mui
            (
                const label speciei,
                const volScalarField& p,
                const volScalarField& T
            ) const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const FluidMulticomponentThermo<BaseThermo>&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

#ifdef NoRepository
    #include "FluidMulticomponentThermo.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
