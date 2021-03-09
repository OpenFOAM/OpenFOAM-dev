/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2021 OpenFOAM Foundation
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
    Foam::fv::semiImplicitSource

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

    Valid fvModels for the \c volumeMode entry include:
    - absolute: values are given as \<quantity\>
    - specific: values are given as \<quantity\>/m3

See also
    Foam::fvModel

SourceFiles
    semiImplicitSource.C

\*---------------------------------------------------------------------------*/

#ifndef semiImplicitSource_H
#define semiImplicitSource_H

#include "fvModel.H"
#include "fvCellSet.H"
#include "HashPtrTable.H"
#include "objectFunction1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
                     Class semiImplicitSource Declaration
\*---------------------------------------------------------------------------*/

class semiImplicitSource
:
    public fvModel
{
public:

    // Public data

        //- Enumeration for volume types
        enum class volumeMode
        {
            absolute,
            specific
        };

        //- Property type names
        static const NamedEnum<volumeMode, 2> volumeModeNames_;


private:

    // Private member data

        //- The set of cells the fvConstraint applies to
        fvCellSet set_;

        //- Volume mode
        volumeMode volumeMode_;

        //- Volume normalisation
        scalar VDash_;

        //- Explicit parts of the sources
        HashPtrTable<objectFunction1> fieldSu_;

        //- Implicit parts of the sources
        HashPtrTable<Function1<scalar>> fieldSp_;


    // Private Member Functions

        //- Non-virtual read
        void readCoeffs();


        // Sources

            //- Add a source term to an equation
            template <class Type>
            void addSupType(fvMatrix<Type>& eqn, const word& fieldName) const;

            //- Add a source term to a compressible equation
            template <class Type>
            void addSupType
            (
                const volScalarField& rho,
                fvMatrix<Type>& eqn,
                const word& fieldName
            ) const;

            //- Add a source term to a phase equation
            template <class Type>
            void addSupType
            (
                const volScalarField& alpha,
                const volScalarField& rho,
                fvMatrix<Type>& eqn,
                const word& fieldName
            ) const;


public:

    //- Runtime type information
    TypeName("semiImplicitSource");


    // Constructors

        //- Construct from components
        semiImplicitSource
        (
            const word& name,
            const word& modelType,
            const dictionary& dict,
            const fvMesh& mesh
        );

        //- Destructor
        virtual ~semiImplicitSource();


    // Member Functions

        // Checks

            //- Return the list of fields for which the fvModel adds source term
            //  to the transport equation
            virtual wordList addSupFields() const;


        // Sources

            //- Add a source term to an equation
            FOR_ALL_FIELD_TYPES(DEFINE_FV_MODEL_ADD_SUP);

            //- Add a source term to a compressible equation
            FOR_ALL_FIELD_TYPES(DEFINE_FV_MODEL_ADD_RHO_SUP);

            //- Add a source term to a phase equation
            FOR_ALL_FIELD_TYPES(DEFINE_FV_MODEL_ADD_ALPHA_RHO_SUP);


        // Mesh motion

            //- Update for mesh changes
            virtual void updateMesh(const mapPolyMesh&);


        // IO

            //- Read source dictionary
            virtual bool read(const dictionary& dict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
