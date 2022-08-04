/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2022 OpenFOAM Foundation
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
    Foam::solvers::multicomponentFluid

Description
    Solver module for steady or transient turbulent flow of compressible
    reacting fluids with optional mesh motion and mesh topology changes.

    Uses the flexible PIMPLE (PISO-SIMPLE) solution for time-resolved and
    pseudo-transient and steady simulations.

    Reference:
    \verbatim
        Greenshields, C. J., & Weller, H. G. (2022).
        Notes on Computational Fluid Dynamics: General Principles.
        CFD Direct Ltd.: Reading, UK.
    \endverbatim

SourceFiles
    multicomponentFluid.C

\*---------------------------------------------------------------------------*/

#ifndef multicomponentFluid_H
#define multicomponentFluid_H

#include "isothermalFluid.H"
#include "fluidMulticomponentThermo.H"
#include "combustionModel.H"
#include "fluidMulticomponentThermophysicalTransportModel.H"
#include "multivariateScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace solvers
{

/*---------------------------------------------------------------------------*\
                          Class multicomponentFluid Declaration
\*---------------------------------------------------------------------------*/

class multicomponentFluid
:
    public isothermalFluid
{

protected:

    // Compositon

        fluidMulticomponentThermo& thermo;

        basicSpecieMixture& composition;
        PtrList<volScalarField>& Y;
        autoPtr<combustionModel> reaction;
        multivariateSurfaceInterpolationScheme<scalar>::fieldTable fields;


    // Thermophysical transport

        autoPtr<fluidMulticomponentThermophysicalTransportModel>
            thermophysicalTransport;


private:

    // Private Member Functions

        //- Set rDeltaT for LTS
        virtual void setRDeltaT();


public:

    //- Runtime type information
    TypeName("multicomponentFluid");


    // Constructors

        //- Construct from region mesh
        multicomponentFluid(fvMesh& mesh);

        //- Disallow default bitwise copy construction
        multicomponentFluid(const multicomponentFluid&) = delete;


    //- Destructor
    virtual ~multicomponentFluid();


    // Member Functions

        //- Construct and solve the energy equation,
        //  convert to temperature
        //  and update thermophysical and transport properties
        virtual void thermophysicalPredictor();

        //- Correct the thermophysical transport modelling
        virtual void thermophysicalTransportCorrector();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const multicomponentFluid&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solvers
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
