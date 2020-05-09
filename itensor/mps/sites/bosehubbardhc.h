//
// Distributed under the ITensor Library License, Version 1.2
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_BOSEHUBBARDHC_H
#define __ITENSOR_BOSEHUBBARDHC_H
#include "itensor/mps/siteset.h"

namespace itensor {

class BoseHubbardHCSite;

using BoseHubbardHC = BasicSiteSet<BoseHubbardHCSite>;

class BoseHubbardHCSite
    {
    IQIndex s;
    public:

    BoseHubbardHCSite() { }

    BoseHubbardHCSite(IQIndex I) : s(I) { }

    BoseHubbardHCSite(int n, Args const& args = Args::global())
        {
        auto conserveNf = args.getBool("ConserveNf",true);
        //auto conserveSz = args.getBool("ConserveSz",true);
        //int Up = (conserveSz ? +1 : 0),
        //    Dn = -Up;
        if(conserveNf)
            {
            s = IQIndex{nameint("site=",n),
                    Index(nameint("0 ",n),1,Site), QN("Nf=",0),
                    Index(nameint("1 ",n),1,Site), QN("Nf=",1)};
            }
        else //don't conserve Nf, only fermion parity
            {
            Error("Nf must be conserved");

           }
        }

    IQIndex
    index() const { return s; }

    IQIndexVal
    state(std::string const& state)
        {
        if(state == "0") 
            {
            return s(1);
            }
        else 
        if(state == "1") 
            {
            return s(2);
            }
        else
            {
            Error("State " + state + " not recognized");
            }
        return IQIndexVal{};
        }

	IQTensor
	op(std::string const& opname,
	   Args const& args) const
        {
        auto sP = prime(s);

        IQIndexVal Zero(s(1)),
                   ZeroP(sP(1)),
                   One(s(2)),
                   OneP(sP(2));

        IQTensor Op(dag(s),sP);

        if(opname == "N")
            {
            Op.set(One,OneP,1);
            }
        else
        if(opname == "NN")
            {
            Op.set(One,OneP,1);
            }
        else
        if(opname == "NInt")
            {
            }
        else
        if(opname == "A")
            {
            Op.set(One,ZeroP,std::sqrt(1.0)); 
            }
        else
        if(opname == "Adag")
            {
            Op.set(Zero,OneP,std::sqrt(1.0)); 
            }
        else
            {
            Error("Operator \"" + opname + "\" name not recognized");
            }

        return Op;
        }
    };


} //namespace itensor

#endif
