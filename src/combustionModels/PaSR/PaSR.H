/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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
    Foam::combustionModels::PaSR

Description
    Partially stirred reactor combustion model.  The model calculates a finite
    rate, based on both turbulence and chemistry time scales.  Depending on
    mesh resolution, the Cmix parameter can be used to scale the turbulence
    mixing time scale.

SourceFiles
    PaSR.C

\*---------------------------------------------------------------------------*/

#ifndef PaSR_H
#define PaSR_H

#include "laminar.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace combustionModels
{

/*---------------------------------------------------------------------------*\
                            Class PaSR Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class PaSR
:
    public laminar<Type>
{
    // Private data

        //- Mixing constant
        scalar Cmix_;

        //- Turbulent reaction switch
        Switch turbulentReaction_;

        //- Mixing parameter
        volScalarField kappa_;


    // Private Member Functions

        //- Disallow copy construct
        PaSR(const PaSR&);

        //- Disallow default bitwise assignment
        void operator=(const PaSR&);


public:

    //- Runtime type information
    TypeName("PaSR");


    // Constructors

        //- Construct from components
        PaSR(const word& modelType, const fvMesh& mesh);


    //- Destructor
    virtual ~PaSR();


    // Member Functions

        // Evolution

            //- Correct combustion rate
            virtual void correct();

            //- Fuel consumption rate matrix.
            virtual tmp<fvScalarMatrix> R(volScalarField& Y) const;

            //- Heat release rate calculated from fuel consumption rate matrix
            virtual tmp<volScalarField> dQ() const;

            //-  Return source for enthalpy equation [kg/m/s3]
            virtual tmp<volScalarField> Sh() const;


    // I-O

            //- Update properties from given dictionary
            virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace combustionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "PaSR.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
