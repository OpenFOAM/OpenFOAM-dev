/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2024 OpenFOAM Foundation
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
    Foam::XiModels::uniformEquilibrium

Description
    Simple uniform equilibrium model for flame wrinkling \c Xi

    See Foam::XiModel for more details on flame wrinkling modelling.

    Reference:
    \verbatim
        Weller, H. G. (1993).
        The development of a new flame area combustion model
        using conditional averaging.
        Thermo-fluids section report TF 9307.
    \endverbatim

SourceFiles
    uniformEquilibrium.C

\*---------------------------------------------------------------------------*/

#ifndef uniformEquilibrium_H
#define uniformEquilibrium_H

#include "XiModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace XiModels
{

/*---------------------------------------------------------------------------*\
                            Class uniformEquilibrium Declaration
\*---------------------------------------------------------------------------*/

class uniformEquilibrium
:
    public XiModel
{
    // Private Data

        //- Equilibrium flame-wrinkling
        scalar XiEq_;


protected:

        //- Update coefficients from given dictionary
        virtual bool readCoeffs(const dictionary& dict);


public:

    //- Runtime type information
    TypeName("uniformEquilibrium");


    // Constructors

        //- Construct from components
        uniformEquilibrium
        (
            const dictionary& dict,
            const psiuMulticomponentThermo& thermo,
            const fluidThermoThermophysicalTransportModel& turbulence,
            const volScalarField& Su
        );

        //- Disallow default bitwise copy construction
        uniformEquilibrium(const uniformEquilibrium&) = delete;


    //- Destructor
    virtual ~uniformEquilibrium();


    // Member Functions

        //- Correct the flame-wrinkling Xi
        virtual void correct()
        {}


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const uniformEquilibrium&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace XiModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
