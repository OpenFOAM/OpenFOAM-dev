/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013 OpenFOAM Foundation
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
    Foam::combustionModels::laminar

Description
    Laminar combustion model.

SourceFiles
    laminar.C

\*---------------------------------------------------------------------------*/

#ifndef laminar_H
#define laminar_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace combustionModels
{

/*---------------------------------------------------------------------------*\
                            Class laminar Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class laminar
:
    public Type
{
    // Private data

        //- Integrate reaction rate over the time-step
        //  using the selected ODE solver
        bool integrateReactionRate_;

protected:

    // Protected Member Functions

        //- Return the chemical time scale
        tmp<volScalarField> tc() const;

private:

    // Private Member Functions

        //- Disallow copy construct
        laminar(const laminar&);

        //- Disallow default bitwise assignment
        void operator=(const laminar&);


public:

    //- Runtime type information
    TypeName("laminar");


    // Constructors

        //- Construct from components
        laminar(const word& modelType, const fvMesh& mesh);


    //- Destructor
    virtual ~laminar();


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
#   include "laminar.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
