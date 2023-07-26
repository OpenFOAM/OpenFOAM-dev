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
    Foam::rhoThermo

Description
    Base-class for thermodynamic properties based on density.

See also
    Foam::basicThermo

SourceFiles
    rhoThermo.C

\*---------------------------------------------------------------------------*/

#ifndef rhoThermo_H
#define rhoThermo_H

#include "basicThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class rhoThermo Declaration
\*---------------------------------------------------------------------------*/

class rhoThermo
:
    virtual public basicThermo
{
public:

    // Public Classes

        //- Forward declare the implementation class
        class implementation;


    //- Destructor
    virtual ~rhoThermo();


    // Member Functions

        // Fields derived from thermodynamic state variables

            //- Density [kg/m^3]
            virtual tmp<volScalarField> rho() const = 0;

            //- Density for patch [kg/m^3]
            virtual tmp<scalarField> rho(const label patchi) const = 0;

            //- Return non-const access to the local density field [kg/m^3]
            virtual volScalarField& rho() = 0;
};


/*---------------------------------------------------------------------------*\
                  Class rhoThermo::implementation Declaration
\*---------------------------------------------------------------------------*/

class rhoThermo::implementation
:
    virtual public rhoThermo
{
protected:

    // Protected data

        //- Density field [kg/m^3]
        volScalarField rho_;


public:

    // Constructors

        //- Construct from dictionary, mesh and phase name
        implementation(const dictionary&, const fvMesh&, const word&);

        //- Disallow default bitwise copy construction
        implementation(const implementation&) = delete;


    //- Destructor
    virtual ~implementation();


    // Member Functions

        // Fields derived from thermodynamic state variables

            //- Density [kg/m^3]
            virtual tmp<volScalarField> rho() const;

            //- Density for patch [kg/m^3]
            virtual tmp<scalarField> rho(const label patchi) const;

            //- Return non-const access to the local density field [kg/m^3]
            virtual volScalarField& rho();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const implementation&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
