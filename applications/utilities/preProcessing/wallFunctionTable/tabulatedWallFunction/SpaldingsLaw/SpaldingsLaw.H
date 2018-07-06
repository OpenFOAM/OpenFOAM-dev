/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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
    Foam::tabulatedWallFunctions::SpaldingsLaw

Description
    Computes U+ as a function of Reynolds number by inverting Spaldings law.

    Example dictionary specification:

        tabulatedWallFunction SpaldingsLaw;

        // Output table info
        tableName       uPlusWallFunctionData; // Output table name
        log10           yes;                // Rey interpreted as log10(Rey)
        dx              0.2;                // Interval log10(Rey)
        x0              -3;                 // Minimum log10(Rey)
        xMax            7;                  // Maximum log10(Rey)

        SpaldingsLawCoeffs
        {
            kappa           0.41;           // Von Karman constant
            E               9.8;            // Law-of-the-wall E coefficient
        }


SourceFiles
    SpaldingsLaw.C

\*---------------------------------------------------------------------------*/

#ifndef SpaldingsLaw_H
#define SpaldingsLaw_H

#include "tabulatedWallFunction.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace tabulatedWallFunctions
{

/*---------------------------------------------------------------------------*\
                        Class SpaldingsLaw Declaration
\*---------------------------------------------------------------------------*/

class SpaldingsLaw
:
    public tabulatedWallFunction
{
protected:

    // Protected data

        //- Von Karman constant
        scalar kappa_;

        //- Law-of-the-wall E coefficient
        scalar E_;


        // Newton iteration solution properties

            //- Maximum number of iterations
            static const label maxIters_;

            //- Tolerance
            static const scalar tolerance_;


    // Protected Member Functions

        //- Invert the function
        virtual void invertFunction();


public:

    //- Run-time type information
    TypeName("SpaldingsLaw");


    // Constructors
    SpaldingsLaw(const dictionary& dict, const polyMesh& mesh);

    //- Destructor
    virtual ~SpaldingsLaw();


    // Member Functions

        // Access

            //- Return y+ as a function of u+
            virtual scalar yPlus(const scalar uPlus) const;

            //- Return Reynolds number as a function of u+
            virtual scalar Re(const scalar uPlus) const;


        // I-O

            //- Write to Ostream
            virtual void writeData(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace tabulatedWallFunctions
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
