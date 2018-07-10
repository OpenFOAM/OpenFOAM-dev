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
    Foam::graph

Description
    Class to create, store and output qgraph files.

SourceFiles
    graph.C

\*---------------------------------------------------------------------------*/

#ifndef graph_H
#define graph_H

#include "string.H"
#include "point.H"
#include "HashPtrTable.H"
#include "curve.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

class graph;

Ostream& operator<<(Ostream&, const graph&);


/*---------------------------------------------------------------------------*\
                           Class graph Declaration
\*---------------------------------------------------------------------------*/

class graph
:
    public HashPtrTable<curve>
{
    // private data

        string title_;
        string xName_;
        string yName_;

        scalarField x_;


        struct xy
        {
            scalar x_, y_;

            xy()
            {}

        // Friend Operators

            friend bool operator==(const xy& a, const xy& b)
            {
                return equal(a.x_, b.x_) && equal(a.y_, b.y_);
            }

            friend bool operator!=(const xy& a, const xy& b)
            {
                return !(a == b);
            }

            friend Istream& operator>>(Istream& is, xy& xyd)
            {
                is >> xyd.x_ >> xyd.y_;
                return is;
            }

            friend Ostream& operator<<(Ostream& os, const xy& xyd)
            {
                os << xyd.x_ << ' ' << xyd.y_;
                return os;
            }
        };


    // Private Member Functions

        void readCurves(Istream&);


public:

    // Constructors

        //- Construct from title and labels (no curves)
        graph
        (
            const string& title,
            const string& xName,
            const string& yName,
            const scalarField& x
        );

        //- Construct from title, labels and y data for 1 curve
        graph
        (
            const string& title,
            const string& xName,
            const string& yName,
            const scalarField& x,
            const scalarField& y
        );

        //- Construct from Istream given title and labels
        graph
        (
            const string& title,
            const string& xName,
            const string& yName,
            Istream& is
        );

        //- Construct from Istream
        graph(Istream& is);


    // Member functions

        // Access

            const string& title() const
            {
                return title_;
            }

            const string& xName() const
            {
                return xName_;
            }

            const string& yName() const
            {
                return yName_;
            }


            const scalarField& x() const
            {
                return x_;
            }

            scalarField& x()
            {
                return x_;
            }


            const scalarField& y() const;

            scalarField& y();


        // Write

            //- Abstract base class for a graph writer
            class writer
            {

            protected:

                void writeXY
                (
                    const scalarField& x,
                    const scalarField& y,
                    Ostream&
                ) const;

            public:

                //- Runtime type information
                TypeName("writer");

                //- Declare run-time constructor selection table
                declareRunTimeSelectionTable
                (
                    autoPtr,
                    writer,
                    word,
                    (),
                    ()
                );


                // Selectors

                    //- Return a reference to the selected writer
                    static autoPtr<writer> New
                    (
                        const word& writeFormat
                    );


                // Constructors

                    //- Construct null
                    writer()
                    {}


                //- Destructor
                virtual ~writer()
                {}


                // Member Functions

                    // Access

                        //- Return the appropriate fileName extension
                        //  for this graph format
                        virtual const word& ext() const = 0;


                    // Write

                        //- Write graph in appropriate format
                        virtual void write(const graph&, Ostream&) const = 0;
            };

            //- Write out graph data as a simple table
            void writeTable(Ostream&) const;

            //- Write graph to stream in given format
            void write(Ostream&, const word& format) const;

            //- Write graph to file in given path-name and format
            void write(const fileName& pName, const word& format) const;

            //- Write graph to file in given path, name and format
            void write
            (
                const fileName& path,
                const word& name,
                const word& format
            ) const;

            //- Helper function to convert string name into appropriate word
            static word wordify(const string& sname);


    // Friend operators

        //- Ostream Operator
        friend Ostream& operator<<(Ostream&, const graph&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
