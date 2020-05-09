//
// Distributed under the ITensor Library License, Version 1.2
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_BOSEHUBBARD5_H
#define __ITENSOR_BOSEHUBBARD5_H
#include "itensor/mps/siteset.h"

namespace itensor {

class BoseHubbard5Site;

using BoseHubbard5 = BasicSiteSet<BoseHubbard5Site>;

class BoseHubbard5Site
    {
    IQIndex s;
    public:

    BoseHubbard5Site() { }

    BoseHubbard5Site(IQIndex I) : s(I) { }

    BoseHubbard5Site(int n, Args const& args = Args::global())
        {
        auto conserveNf = args.getBool("ConserveNf",true);
        //auto conserveSz = args.getBool("ConserveSz",true);
        //int Up = (conserveSz ? +1 : 0),
        //    Dn = -Up;
        if(conserveNf)
            {
            s = IQIndex{nameint("site=",n),
                    Index(nameint("0 ",n),1,Site), QN("Nf=",0),
                    Index(nameint("1 ",n),1,Site), QN("Nf=",1),
                    Index(nameint("2 ",n),1,Site), QN("Nf=",2),
                    Index(nameint("3 ",n),1,Site), QN("Nf=",3),
                    Index(nameint("4 ",n),1,Site), QN("Nf=",4),
                    Index(nameint("5 ",n),1,Site), QN("Nf=",5)};
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
        if(state == "2") 
            {
            return s(3);
            }
        else 
        if(state == "3") 
            {
            return s(4);
            }
        else 
        if(state == "4") 
            {
            return s(5);
            }
        else 
        if(state == "5") 
            {
            return s(6);
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
                   OneP(sP(2)),
                   Two(s(3)),
                   TwoP(sP(3)),
                   Three(s(4)),
                   ThreeP(sP(4)),
                   Four(s(5)),
                   FourP(sP(5)),
                   Five(s(6)),
                   FiveP(sP(6));

        IQTensor Op(dag(s),sP);

        if(opname == "N")
            {
            Op.set(One,OneP,1);
            Op.set(Two,TwoP,2);
            Op.set(Three,ThreeP,3);
            Op.set(Four,FourP,4);
            Op.set(Five,FiveP,5);
            }
        else
        if(opname == "NN")
            {
            Op.set(One,OneP,1);
            Op.set(Two,TwoP,4);
            Op.set(Three,ThreeP,9);
            Op.set(Four,FourP,16);
            Op.set(Five,FiveP,25);
            }
        else
        if(opname == "NInt")
            {
            Op.set(Two,TwoP,2);
            Op.set(Three,ThreeP,6);
            Op.set(Four,FourP,12);
            Op.set(Five,FiveP,20);
            }
        else
        if(opname == "A")
            {
            Op.set(One,ZeroP,std::sqrt(1.0)); 
            Op.set(Two,OneP,std::sqrt(2.0)); 
            Op.set(Three,TwoP,std::sqrt(3.0));
            Op.set(Four,ThreeP,std::sqrt(4.0));
            Op.set(Five,FourP,std::sqrt(5.0));
            }
        else
        if(opname == "Adag")
            {
            Op.set(Zero,OneP,std::sqrt(1.0)); 
            Op.set(One,TwoP,std::sqrt(2.0));
            Op.set(Two,ThreeP,std::sqrt(3.0));
            Op.set(Three,FourP,std::sqrt(4.0));
            Op.set(Four,FiveP,std::sqrt(5.0));
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
