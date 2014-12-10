
template<int N, int I>
class VectorSpaceOps
{
public:

    static const int endLoop = (I < N-1) ? 1 : 0;

    template<class V, class S, class EqOp>
    static inline void eqOpS(V& vs, const S& s, EqOp eo)
    {
        eo(vs.v_[I], s);
        VectorSpaceOps<endLoop*N, endLoop*(I+1)>::eqOpS(vs, s, eo);
    }

    template<class S, class V, class EqOp>
    static inline void SeqOp(S& s, const V& vs, EqOp eo)
    {
        eo(s, vs.v_[I]);
        VectorSpaceOps<endLoop*N, endLoop*(I+1)>::SeqOp(s, vs, eo);
    }

    template<class V1, class V2, class EqOp>
    static inline void eqOp(V1& vs1, const V2& vs2, EqOp eo)
    {
        eo(vs1.v_[I], vs2.v_[I]);
        VectorSpaceOps<endLoop*N, endLoop*(I+1)>::eqOp(vs1, vs2, eo);
    }


    template<class V, class V1, class S, class Op>
    static inline void opVS(V& vs, const V1& vs1, const S& s, Op o)
    {
        vs.v_[I] = o(vs1.v_[I], s);
        VectorSpaceOps<endLoop*N, endLoop*(I+1)>::opVS(vs, vs1, s, o);
    }

    template<class V, class S, class V1, class Op>
    static inline void opSV(V& vs, const S& s,  const V1& vs1, Op o)
    {
        vs.v_[I] = o(s, vs1.v_[I]);
        VectorSpaceOps<endLoop*N, endLoop*(I+1)>::opSV(vs, s, vs1, o);
    }

    template<class V, class V1, class Op>
    static inline void op(V& vs, const V1& vs1, const V1& vs2, Op o)
    {
        vs.v_[I] = o(vs1.v_[I], vs2.v_[I]);
        VectorSpaceOps<endLoop*N, endLoop*(I+1)>::op(vs, vs1, vs2, o);
    }
};


template<>
class VectorSpaceOps<0, 0>
{
public:

    template<class V, class S, class EqOp>
    static inline void eqOpS(V&, const S&, EqOp)
    {}

    template<class S, class V, class EqOp>
    static inline void SeqOp(S&, const V&, EqOp)
    {}

    template<class V1, class V2, class EqOp>
    static inline void eqOp(V1&, const V2&, EqOp)
    {}


    template<class V, class V1, class S, class Op>
    static inline void opVS(V& vs, const V1&, const S&, Op)
    {}

    template<class V, class S, class V1, class Op>
    static inline void opSV(V& vs, const S&, const V1&, Op)
    {}

    template<class V, class V1, class Op>
    static inline void op(V& vs, const V1&, const V1&, Op)
    {}
};
