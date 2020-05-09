//
// Distributed under the ITensor Library License, Version 1.2
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_BOSEHUBBARDGEN_H
#define __ITENSOR_BOSEHUBBARDGEN_H
#include "itensor/mps/siteset.h"

namespace itensor {

class BoseHubbardGenSite;

using BoseHubbardGen = BasicSiteSet<BoseHubbardGenSite>;

class BoseHubbardGenSite
    {
    IQIndex s;


    public:
    static int minnp;
    static int rangenp;


    BoseHubbardGenSite() { }

    BoseHubbardGenSite(IQIndex I) : s(I) { }

    BoseHubbardGenSite(int n, Args const& args = Args::global())
        {
        auto conserveNf = args.getBool("ConserveNf",true);
	minnp=args.getInt("minnp",0);
	rangenp=args.getInt("rangenp",4);
        //auto conserveSz = args.getBool("ConserveSz",true);
        //int Up = (conserveSz ? +1 : 0),
        //    Dn = -Up;
        if(conserveNf)
            {
	    if(rangenp==1){
		    s = IQIndex{nameint("site=",n),
		            Index(nameint("0 ",n),1,Site), QN("Nf=",0),
		            Index(nameint("1 ",n),1,Site), QN("Nf=",1)};
	    }else
	    if(rangenp==2){
		    s = IQIndex{nameint("site=",n),
		            Index(nameint("0 ",n),1,Site), QN("Nf=",0),
		            Index(nameint("1 ",n),1,Site), QN("Nf=",1),
		            Index(nameint("2 ",n),1,Site), QN("Nf=",2)};
	    }else
	    if(rangenp==3){
		    s = IQIndex{nameint("site=",n),
		            Index(nameint("0 ",n),1,Site), QN("Nf=",0),
		            Index(nameint("1 ",n),1,Site), QN("Nf=",1),
		            Index(nameint("2 ",n),1,Site), QN("Nf=",2),
		            Index(nameint("3 ",n),1,Site), QN("Nf=",3)};
	    }else
	    if(rangenp==4){
		    s = IQIndex{nameint("site=",n),
		            Index(nameint("0 ",n),1,Site), QN("Nf=",0),
		            Index(nameint("1 ",n),1,Site), QN("Nf=",1),
		            Index(nameint("2 ",n),1,Site), QN("Nf=",2),
		            Index(nameint("3 ",n),1,Site), QN("Nf=",3),
		            Index(nameint("4 ",n),1,Site), QN("Nf=",4)};
	    }else
	    if(rangenp==5){
		    s = IQIndex{nameint("site=",n),
		            Index(nameint("0 ",n),1,Site), QN("Nf=",0),
		            Index(nameint("1 ",n),1,Site), QN("Nf=",1),
		            Index(nameint("2 ",n),1,Site), QN("Nf=",2),
		            Index(nameint("3 ",n),1,Site), QN("Nf=",3),
		            Index(nameint("4 ",n),1,Site), QN("Nf=",4),
			    Index(nameint("5 ",n),1,Site), QN("Nf=",5)};
	    }else
	    if(rangenp==6){
		    s = IQIndex{nameint("site=",n),
		            Index(nameint("0 ",n),1,Site), QN("Nf=",0),
		            Index(nameint("1 ",n),1,Site), QN("Nf=",1),
		            Index(nameint("2 ",n),1,Site), QN("Nf=",2),
		            Index(nameint("3 ",n),1,Site), QN("Nf=",3),
		            Index(nameint("4 ",n),1,Site), QN("Nf=",4),
			    Index(nameint("5 ",n),1,Site), QN("Nf=",5),
			    Index(nameint("6 ",n),1,Site), QN("Nf=",6)};
	    }else
	    if(rangenp==7){
		    s = IQIndex{nameint("site=",n),
		            Index(nameint("0 ",n),1,Site), QN("Nf=",0),
		            Index(nameint("1 ",n),1,Site), QN("Nf=",1),
		            Index(nameint("2 ",n),1,Site), QN("Nf=",2),
		            Index(nameint("3 ",n),1,Site), QN("Nf=",3),
		            Index(nameint("4 ",n),1,Site), QN("Nf=",4),
			    Index(nameint("5 ",n),1,Site), QN("Nf=",5),
			    Index(nameint("6 ",n),1,Site), QN("Nf=",6),
			    Index(nameint("7 ",n),1,Site), QN("Nf=",7)};
	    }else
	    if(rangenp==8){
		    s = IQIndex{nameint("site=",n),
		            Index(nameint("0 ",n),1,Site), QN("Nf=",0),
		            Index(nameint("1 ",n),1,Site), QN("Nf=",1),
		            Index(nameint("2 ",n),1,Site), QN("Nf=",2),
		            Index(nameint("3 ",n),1,Site), QN("Nf=",3),
		            Index(nameint("4 ",n),1,Site), QN("Nf=",4),
			    Index(nameint("5 ",n),1,Site), QN("Nf=",5),
			    Index(nameint("6 ",n),1,Site), QN("Nf=",6),
			    Index(nameint("7 ",n),1,Site), QN("Nf=",7),
			    Index(nameint("8 ",n),1,Site), QN("Nf=",8)};
	    }else
	    if(rangenp==9){
		    s = IQIndex{nameint("site=",n),
		            Index(nameint("0 ",n),1,Site), QN("Nf=",0),
		            Index(nameint("1 ",n),1,Site), QN("Nf=",1),
		            Index(nameint("2 ",n),1,Site), QN("Nf=",2),
		            Index(nameint("3 ",n),1,Site), QN("Nf=",3),
		            Index(nameint("4 ",n),1,Site), QN("Nf=",4),
			    Index(nameint("5 ",n),1,Site), QN("Nf=",5),
			    Index(nameint("6 ",n),1,Site), QN("Nf=",6),
			    Index(nameint("7 ",n),1,Site), QN("Nf=",7),
			    Index(nameint("8 ",n),1,Site), QN("Nf=",8),
			    Index(nameint("9 ",n),1,Site), QN("Nf=",9)};
	    }else
	    if(rangenp==10){
		    s = IQIndex{nameint("site=",n),
		            Index(nameint("0 ",n),1,Site), QN("Nf=",0),
		            Index(nameint("1 ",n),1,Site), QN("Nf=",1),
		            Index(nameint("2 ",n),1,Site), QN("Nf=",2),
		            Index(nameint("3 ",n),1,Site), QN("Nf=",3),
		            Index(nameint("4 ",n),1,Site), QN("Nf=",4),
			    Index(nameint("5 ",n),1,Site), QN("Nf=",5),
			    Index(nameint("6 ",n),1,Site), QN("Nf=",6),
			    Index(nameint("7 ",n),1,Site), QN("Nf=",7),
			    Index(nameint("8 ",n),1,Site), QN("Nf=",8),
			    Index(nameint("9 ",n),1,Site), QN("Nf=",9),
			    Index(nameint("10 ",n),1,Site), QN("Nf=",10)};
	    }else{
            	Error("rangenp "+ std::to_string(rangenp)+" is not supported");
	    }

        }else //don't conserve Nf, only fermion parity
            {
            Error("Nf must be conserved");

           }
        }

    IQIndex
    index() const { return s; }

    IQIndexVal
    state(std::string const& state)
        {
	int statenum=std::stoi(state);
	if(statenum>rangenp){
            Error("State "+std::to_string(statenum)+" larger than "+std::to_string(rangenp) );
            return IQIndexVal{};
	}else{
	   return s(statenum+1);
	}
	

        }

	IQTensor
	op(std::string const& opname,
	   Args const& args) const
        {
        auto sP = prime(s);

        //println(rangenp," ",minnp);
        IQTensor Op(dag(s),sP);

        if(opname == "N")
            {
	    for(int b = 0; b <= rangenp; ++b){
		Op.set(s(b+1),sP(b+1),b+minnp);
	    }

            }
        else
        if(opname == "NN")
            {
	    for(int b = 0; b <= rangenp; ++b){
		Op.set(s(b+1),sP(b+1),(b+minnp)*(b+minnp));
	    }
            }
        else
        if(opname == "NInt")
            {
	    for(int b = 0; b <= rangenp; ++b){
		Op.set(s(b+1),sP(b+1),(b+minnp)*(b+minnp-1));
	    }
            }
        else
        if(opname == "A")
            {
	    for(int b = 1; b <= rangenp; ++b){
		Op.set(s(b+1),sP(b-1+1),std::sqrt(b+minnp));
	    }
            }
        else
        if(opname == "Adag")
            {
	    for(int b = 0; b < rangenp; ++b){
		Op.set(s(b+1),sP(b+1+1),std::sqrt(b+1+minnp));
	    }
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
