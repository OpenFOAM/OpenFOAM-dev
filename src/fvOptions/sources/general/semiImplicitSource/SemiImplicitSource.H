/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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
    Foam::fv::SemiImplicitSource

Description
    Semi-implicit source, described using an input dictionary.  The injection
    rate coefficients are specified as pairs of Su-Sp coefficients, i.e.

        \f[
            S(x) = S_u + S_p x
        \f]

    where
    \vartable
        S(x)    | net source for field 'x'
        S_u     | explicit source contribution
        S_p     | linearised implicit contribution
    \endvartable

    Example tabulated heat source specification for internal energy:
    \verbatim
    volumeMode      absolute; // specific
    sources
    {
        e
        {
            explicit table ((0 0) (1.5 $power));
            implicit 0;
        }
    }
    \endverbatim

    Example coded heat source specification for enthalpy:
    \verbatim
    volumeMode      absolute; // specific
    sources
    {
        h
        {
            explicit
            {
                type coded;
                name heatInjection;
                code
                #{
                    // Power amplitude
                    const scalar powerAmplitude = 1000;

                    // x is the current time
                    return mag(powerAmplitude*sin(x));
                #};
            }
            implicit 0;
        }
    }
    \endverbatim

    Valid options for the \c volumeMode entry include:
    - absolute: values are given as \<quantity\>
    - specific: values are given as \<quantity\>/m3

See also
    Foam::fvOption

SourceFiles
    SemiImplicitSource.C

\*---------------------------------------------------------------------------*/

#ifndef SemiImplicitSource_H
#define SemiImplicitSource_H

#include "cellSetOption.H"
#include "Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
                     Class SemiImplicitSource Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class SemiImplicitSource
:
    public cellSetOption
{
public:

    // Public data

        //- Enumeration for volume types
        enum volumeModeType
        {
            vmAbsolute,
            vmSpecific
        };


private:

    // Private member data

        //- Word list of volume mode type names
        static const wordList volumeModeTypeNames_;

        //- Volume mode
        volumeModeType volumeMode_;

        //- Volume normalisation
        scalar VDash_;

        //- Explicit source functions for the fields
        PtrList<Function1<Type>> Su_;

        //- Implicit source functions for the fields
        PtrList<Function1<scalar>> Sp_;


    // Private member functions

        //- Helper function to convert from a word to a volumeModeType
        volumeModeType wordToVolumeModeType(const word& vtName) const;

        //- Helper function to convert from a volumeModeType to a word
        word volumeModeTypeToWord(const volumeModeType& vtType) const;

        //- Set the local field data
        void setFieldData(const dictionary& dict);


public:

    //- Runtime type information
    TypeName("SemiImplicitSource");


    // Constructors

        //- Construct from components
        SemiImplicitSource
        (
            const word& name,
            const word& modelType,
            const dictionary& dict,
            const fvMesh& mesh
        );


    // Member Functions

        // Evaluation

            //- Add explicit contribution to equation
            virtual void addSup
            (
                fvMatrix<Type>& eqn,
                const label fieldi
            );

            //- Add explicit contribution to compressible equation
            virtual void addSup
            (
                const volScalarField& rho,
                fvMatrix<Type>& eqn,
                const label fieldi
            );


        // IO

            //- Read source dictionary
            virtual bool read(const dictionary& dict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "SemiImplicitSource.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
