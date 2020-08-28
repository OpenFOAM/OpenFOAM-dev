/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2017-2020 OpenFOAM Foundation
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

Description:
    Functions which cut triangles and tetrahedra. Generic operations are
    applied to each half of a cut.

SourceFiles:
    cutI.H
    cutTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef cut_H
#define cut_H

#include "FixedList.H"
#include "nil.H"
#include "plane.H"
#include "tetPointRef.H"
#include "triPointRef.H"
#include "zero.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace cut
{

/*---------------------------------------------------------------------------*\
                          Class uniformOp Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class uniformOp
{
    //- Data
    Type data_;


public:

    // Constructors

        //- Construct null
        uniformOp()
        {}

        //- Construct from data
        uniformOp(Type data)
        :
            data_(data)
        {}


    // Member Functions

        //- Access the data
        Type data() const
        {
            return data_;
        }
};


/*---------------------------------------------------------------------------*\
                            Class noOp Declaration
\*---------------------------------------------------------------------------*/

class noOp
:
    public uniformOp<nil>
{
public:

    // Typedefs

        //- Result type
        typedef zero result;


    // Constructors

        //- Construct null
        noOp()
        {}

        //- Construct from base
        noOp(const uniformOp<nil>& op)
        :
            uniformOp<nil>(op)
        {}


    // Member Operators

        //- Operate on nothing
        result operator()() const
        {
            return result();
        }

        //- Operate on a triangle or tetrahedron
        template<class Type, unsigned Size>
        result operator()(const FixedList<Type, Size>& p) const
        {
            return result();
        }
};


/*---------------------------------------------------------------------------*\
                           Class areaOp Declaration
\*---------------------------------------------------------------------------*/

class areaOp
:
    public uniformOp<nil>
{
public:

    // Typedefs

        //- Result type
        typedef vector result;


    // Constructors

        //- Construct null
        areaOp()
        {}

        //- Construct from base
        areaOp(const uniformOp<nil>& op)
        :
            uniformOp<nil>(op)
        {}


    // Member Operators

        //- Operate on nothing
        result operator()() const
        {
            return vector::zero;
        }

        //- Operate on a triangle
        result operator()(const FixedList<point, 3>& p) const
        {
            return result(triPointRef(p[0], p[1], p[2]).area());
        }
};


/*---------------------------------------------------------------------------*\
                          Class areaMagOp Declaration
\*---------------------------------------------------------------------------*/

class areaMagOp
:
    public uniformOp<nil>
{
public:

    // Typedefs

        //- Result type
        typedef scalar result;


    // Constructors

        //- Construct null
        areaMagOp()
        {}

        //- Construct from base
        areaMagOp(const uniformOp<nil>& op)
        :
            uniformOp<nil>(op)
        {}


    // Member Operators

        //- Operate on nothing
        result operator()() const
        {
            return 0;
        }

        //- Operate on a triangle
        result operator()(const FixedList<point, 3>& p) const
        {
            return result(triPointRef(p[0], p[1], p[2]).mag());
        }
};


/*---------------------------------------------------------------------------*\
                          Class volumeOp Declaration
\*---------------------------------------------------------------------------*/

class volumeOp
:
    public uniformOp<nil>
{
public:

    // Typedefs

        //- Result type
        typedef scalar result;


    // Constructors

        //- Construct null
        volumeOp()
        {}

        //- Construct from base
        volumeOp(const uniformOp<nil>& op)
        :
            uniformOp<nil>(op)
        {}


    // Member Operators

        //- Operate on nothing
        result operator()() const
        {
            return 0;
        }

        //- Operate on a tetrahedron
        result operator()(const FixedList<point, 4>& p) const
        {
            return result(tetPointRef(p[0], p[1], p[2], p[3]).mag());
        }
};


/*---------------------------------------------------------------------------*\
                       Class areaIntegrateOp Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class areaIntegrateOp
:
    public FixedList<Type, 3>
{
public:

    // Typedefs

        //- Result type
        typedef typename outerProduct<Type, vector>::type result;


    // Constructors

        //- Construct from base
        areaIntegrateOp(const FixedList<Type, 3>& x)
        :
            FixedList<Type, 3>(x)
        {}


    // Member Operators

        //- Operate on nothing
        result operator()() const
        {
            return pTraits<result>::zero;
        }

        //- Operate on a triangle
        result operator()(const FixedList<point, 3>& p) const
        {
            const FixedList<Type, 3>& x = *this;
            return result(areaOp()(p)*(x[0] + x[1] + x[2])/3);
        }
};


/*---------------------------------------------------------------------------*\
                     Class areaMagIntegrateOp Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class areaMagIntegrateOp
:
    public FixedList<Type, 3>
{
public:

    // Typedefs

        //- Result type
        typedef Type result;


    // Constructors

        //- Construct from base
        areaMagIntegrateOp(const FixedList<Type, 3>& x)
        :
            FixedList<Type, 3>(x)
        {}


    // Member Operators

        //- Operate on nothing
        result operator()() const
        {
            return pTraits<result>::zero;
        }

        //- Operate on a triangle
        result operator()(const FixedList<point, 3>& p) const
        {
            const FixedList<Type, 3>& x = *this;
            return result(areaMagOp()(p)*(x[0] + x[1] + x[2])/3);
        }
};


/*---------------------------------------------------------------------------*\
                      Class volumeIntegrateOp Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class volumeIntegrateOp
:
    public FixedList<Type, 4>
{
public:

    // Typedefs

        //- Result type
        typedef Type result;


    // Constructors

        //- Construct from base
        volumeIntegrateOp(const FixedList<Type, 4>& x)
        :
            FixedList<Type, 4>(x)
        {}


    // Member Operators

        //- Operate on nothing
        result operator()() const
        {
            return pTraits<result>::zero;
        }

        //- Operate on a tetrahedron
        result operator()(const FixedList<point, 4>& p) const
        {
            const FixedList<Type, 4>& x = *this;
            return result(volumeOp()(p)*(x[0] + x[1] + x[2] + x[3])/4);
        }
};


/*---------------------------------------------------------------------------*\
                           Class listOp Declaration
\*---------------------------------------------------------------------------*/

template<class Type, unsigned Size>
class listOp
:
    public uniformOp<nil>
{
public:

    // Classes

        //- Result class
        class result
        :
            public DynamicList<FixedList<Type, Size>>
        {
        public:

            // Constructors

                //- Construct from a single element
                result(const FixedList<Type, Size>& x)
                :
                    DynamicList<FixedList<Type, Size>>(1, x)
                {}


            // Member Operators

                //- Add together two lists
                result operator+(const result& x) const
                {
                    result r(*this);
                    r.append(x);
                    return r;
                }
        };


    // Constructors

        //- Construct null
        listOp()
        {}

        //- Construct from base
        listOp(const uniformOp<nil>& op)
        :
            uniformOp<nil>(op)
        {}


    // Member Operators

        //- Operate on nothing
        result operator()() const
        {
            return result();
        }

        //- Operate on a triangle or tetrahedron
        result operator()(const FixedList<Type, Size>& p) const
        {
            return result(p);
        }
};


typedef listOp<point, 3> listTriOp;


typedef listOp<point, 4> listTetOp;


/*---------------------------------------------------------------------------*\
                          Class appendOp Declaration
\*---------------------------------------------------------------------------*/

template<class Container>
class appendOp
:
    public uniformOp<Container&>
{
public:

    // Typedefs

        //- Result type
        typedef zero result;


    // Constructors

        //- Construct from a container reference
        appendOp(Container& x)
        :
            uniformOp<Container&>(x)
        {}

        //- Construct from base
        appendOp(const uniformOp<Container&>& op)
        :
            uniformOp<Container&>(op)
        {}


    // Member Operators

        //- Operate on nothing
        result operator()() const
        {
            return result();
        }

        //- Operate on a triangle or tetrahedron
        template<class Type, unsigned Size>
        result operator()(const FixedList<Type, Size>& p) const
        {
            this->data().append(p);
            return result();
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- Trait to determine the result of the addition of two operations
template<class AheadOp, class BehindOp>
class opAddResult;

template<class Op>
class opAddResult<Op, Op>
{
public:

    typedef typename Op::result type;
};

template<>
class opAddResult<noOp, noOp>
{
public:

    typedef typename noOp::result type;
};

template<class Op>
class opAddResult<noOp, Op>
{
public:

    typedef typename Op::result type;
};

template<class Op>
class opAddResult<Op, noOp>
{
public:

    typedef typename Op::result type;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace cut

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- Cut a triangle along the zero plane defined by the given levels. Templated
//  on aboveOp and belowOp; the operations applied to the regions on either side
//  of the cut.
template<class Point, class AboveOp, class BelowOp>
typename cut::opAddResult<AboveOp, BelowOp>::type triCut
(
    const FixedList<Point, 3>& tri,
    const FixedList<scalar, 3>& level,
    const AboveOp& aboveOp,
    const BelowOp& belowOp
);

//- As above, but with a plane specifying the location of the cut
template<class AboveOp, class BelowOp>
typename cut::opAddResult<AboveOp, BelowOp>::type triCut
(
    const FixedList<point, 3>& tri,
    const plane& s,
    const AboveOp& aboveOp,
    const BelowOp& belowOp
);

//- As triCut, but for a tetrahedron.
template<class Point, class AboveOp, class BelowOp>
typename cut::opAddResult<AboveOp, BelowOp>::type tetCut
(
    const FixedList<Point, 4>& tet,
    const FixedList<scalar, 4>& level,
    const AboveOp& aboveOp,
    const BelowOp& belowOp
);

//- As above, but with a plane specifying the location of the cut
template<class AboveOp, class BelowOp>
typename cut::opAddResult<AboveOp, BelowOp>::type tetCut
(
    const FixedList<point, 4>& tet,
    const plane& s,
    const AboveOp& aboveOp,
    const BelowOp& belowOp
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "cutI.H"

#ifdef NoRepository
    #include "cutTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
