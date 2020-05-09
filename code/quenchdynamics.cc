#include "itensor/all.h"
#include <stdio.h>
#include <complex.h>
#include <iostream>
#include <fstream>
#include <ctime>
#include <algorithm>


/*
Code to simulate quench dynamics in Ring-lead and Y-junctions. Supports bosons and fermions.
First, it calculates ground state via DMRG. You can add a density perturbation in the source of the leads using the potential variable.
Then, it calculates the dynamics. The potential is quenched at t=0 and the resulting dynamics is simulated. The density and current in time is outputed to file.
of output txt: 
Line1: time where snapshots where taken. In each line of length numberSnaps
Line2: current bond dimension
Line3: computational time
Line4: current 1 measured
Line5: current 2 measured
Line6: current 3 measured
Line7: Overlap with initial state
Line8: Not used
Line9 to Number of sites +9: density for each site

*/


using namespace itensor;
using namespace std;

//Efficient implementation of ring Hamiltonian

int BoseHubbardGenSite::minnp;
int BoseHubbardGenSite::rangenp;



class hblueprint //Create hamiltonian from local potential and non-local hopping terms
{
	struct nonlocalterms {
	  int start;
	  int end;
	  string opstart;
	  string opend;
	  complex<double> value;
	};

	struct localterms {
	  int pos;
	  string op;
	  double value;
	};

    public:
    vector<localterms> local;
    vector<nonlocalterms> nonlocal;
     
    //Default Constructor
    hblueprint()
    {

    }
     
    //Parametrized Constructor
    void add(double value, string op, int pos)
    {
	localterms temp;
	temp.pos=pos;
	temp.op=op;
	temp.value=value;
	local.push_back(temp);
    }

    void add(complex<double> value, string opstart,int start, string opend, int end, bool conjugate=false)
    {
	/*if(opstart>opend){ 
		throw std::invalid_argument("startpoint must be smaller than endpoint");
	}*/
	nonlocalterms temp;
	temp.start=start;
	temp.end=end;
	temp.opstart=opstart;
	temp.opend=opend;
	temp.value=value;

	nonlocal.push_back(temp);
	if(conjugate==true){
		nonlocalterms tempconj;
		tempconj.start=end;
		tempconj.end=start;
		tempconj.opstart=opstart;
		tempconj.opend=opend;
		tempconj.value=conj(value);
		nonlocal.push_back(tempconj);
	}
    }

    void print(){
	println("Local terms");
	for (auto it = local.begin(); it!=local.end(); ++it) {
		println("",it->value,", ",it->pos, ", ",it->op);
	}
	println("Non-local term");
	for (auto it = nonlocal.begin(); it!=nonlocal.end(); ++it) {
		println("",it->value,", ",it->start, ", ",it->end,", ",it->opstart, ", ",it->opend);	
	}
    }
	


    IQMPO generateMPO(SiteSet const& sites,complex<double> tau=complex<double> (0.0,0.0))
    {
	auto ampo = AutoMPO(sites);
	for (auto it = local.begin(); it!=local.end(); ++it) {
		if(it->value!=0) ampo += it->value,it->op,it->pos;
	}
	
	for (auto it = nonlocal.begin(); it!=nonlocal.end(); ++it) {
		if(it->value!=complex<double> (0.0,0.0))
		{
			if(it->opstart=="C") throw std::invalid_argument("Start is operator C, this may not be properly implemented (because of fermionic rules, a minus may have to be included");
			ampo += it->value,it->opstart,it->start,it->opend,it->end;

		}
		
	}
	if(tau==complex<double> (0.0,0.0)) return IQMPO(ampo); //make MPO if no tau given (e.g. for expectation values)
	return toExpH<IQTensor>(ampo,tau); //Return MPO for evolving by time tau
    }




};



//Load state from file
void loadsites(SiteSet& sites,string file, int mpstype, Args const& args = Args::global()){
	if(mpstype==0){
		sites=readFromFile<BoseHubbard>(file);
	}else if(mpstype==1 or mpstype==2){
		sites=readFromFile<Spinless>(file);
	}else if(mpstype==3 or mpstype==6){
		sites=readFromFile<Hubbard>(file);
	}else if(mpstype==4){
		sites=readFromFile<BoseHubbard3>(file);
	}else if(mpstype==5){
		sites=readFromFile<BoseHubbard5>(file);
	}else if(mpstype==7){ //Bose Hubbard with arbitrary particles per site, use variable below to set.
		sites=readFromFile<BoseHubbardGen>(file);
		BoseHubbardGenSite::minnp = args.getInt("minnp",0);
		BoseHubbardGenSite::rangenp = args.getInt("rangenp",4);
		//BasicSiteSet<itensor::BoseHubbardGenSite> sitesB=readFromFile<BoseHubbardGen>(file);
		//sitesB::minnp=args.getInt("minnp",0);
		//BoseHubbardGenSite sitesB=(BoseHubbardGenSite)sites;
		//sitesB.minnp=args.getInt("minnp",0);
		//sitesB.rangenp=args.getInt("rangenp",4);
	
	}
}



string ftos(float f, int nd) {
   ostringstream ostr;
   int tens = pow(10,nd);
   ostr << round(f*tens)/tens;
   return ostr.str();
}



string ReplaceAll(string str, const std::string& from, const std::string& to) {
    size_t start_pos = 0;
    while((start_pos = str.find(from, start_pos)) != std::string::npos) {
        str.replace(start_pos, from.length(), to);
        start_pos += to.length(); // Handles case where 'to' is a substring of 'from'
    }
    return str;
}

int cleanupoutput(const std::string& dataset, const std::string& prefix){
   int result=0;
   result+=remove(("psi"+dataset).c_str());
   result+=remove(("sites"+dataset).c_str());
   result+=remove(("Density"+dataset+prefix+".dat").c_str());
   if(result!=0) println("Could not clean up files for " +dataset);

   return result;
}

//Get local expectation value
double getlocalExpectation(int a, const std::string& op,SiteSet const& sites, IQMPS psi){
	return (dag(prime(psi.A(a),Site))*sites.op(op,a)*psi.A(a)).real();
}


//Get expectation value of two diferent states , e.g. <psia*op*psib>
complex<double> getlocalExpectationOverlap(int a, const std::string& op,SiteSet const& sites, IQMPS psia,IQMPS psib){
	//psia.position(a);
	//psib.position(a);
	//return (dag(prime(psib.A(a),Site))*sites.op(op,a)*psia.A(a)).cplx();
	//println("Local",a);
	auto ampo = AutoMPO(sites);
	ampo += 1,op,a;
        std::cout.setstate(std::ios_base::failbit);//Supress spam messages
	auto opMPO = IQMPO(ampo);
        std::cout.clear();
	//println("Local end",a);
	auto result= overlapC(psia,opMPO,psib);
	return result;
}


//Calculates Adag(a)*a(b) for two different psi. For fermions, take care to use C operators, to take care of Jordan Wigner strings
std::complex<double> getcorrelatorOverlap(int a,int b,const std::string& createOp, const std::string& destroyOp,SiteSet const& sites,IQMPS psia,IQMPS psib) {
	//println("Start",a," ",b);
	auto ampo = AutoMPO(sites);
	ampo+=1,createOp,a,destroyOp,b;
        std::cout.setstate(std::ios_base::failbit);
	auto opMPO = IQMPO(ampo);
        std::cout.clear();
	//println("end",a," ",b);
	return overlapC(psia,opMPO,psib);
}


//Calculates Adag(a)*a(b), only for bosons and fermions when nearest-neighbors so far!! Fermion with more than nearest neighbor need jordan-wigner string
std::complex<double> getcorrelator(int a,int b,const std::string& createOp, const std::string& destroyOp,SiteSet const& sites,IQMPS psi) {
	//Given an MPS or IQMPS called "psi",
	//constructed from a SiteSet "sites"

	//Replace "Op1" and "Op2" with the actual names
	//of the operators you want to measure

	auto i=a;
	auto j=b;
	auto op_i = sites.op(createOp,i); //This is for i<j
	auto op_j = sites.op(destroyOp,j);


	if(j<i){
		i=b;
		j=a;
		op_i = sites.op(destroyOp,i);
		op_j = sites.op(createOp,j);
	}

	//below we will assume j > i

	//'gauge' the MPS to site i
	//any 'position' between i and j, inclusive, would work here
	psi.position(i); 

	//psi.Anc(1) *= psi.A(0); //Uncomment if doing iDMRG calculation

	//index linking i to i+1:
	auto ir = commonIndex(psi.A(i),psi.A(i+1),Link);

	auto C = psi.A(i)*op_i*dag(prime(psi.A(i),Site,ir));
	for(int k = i+1; k < j; ++k)
		{
		C *= psi.A(k);
		C *= dag(prime(psi.A(k),Link));
	   

		}
	C *= psi.A(j);
	C *= op_j;
	//index linking j to j-1:
	auto jl = commonIndex(psi.A(j),psi.A(j-1),Link);
	C *= dag(prime(psi.A(j),jl,Site));

	auto result = C.cplx(); //or C.cplx() if expecting complex
	return result;
}


int main(int argc, char* argv[])
    {
    std::clock_t start;
    double duration;


    start = std::clock();

    //auto mpstype=input.getInt("mpstype",0); //which type of mps, 0:Bose-Hubbard with 4 particles per site, 1: hard-core BH, 2: spinless fermions, 3: Spinfull fermions, 4: Bose-Hubbard with 3 particles per site, 5: Bose-Hubbard with 5 particles per site, 6: Hard-core BH two species
    //#define MPSTYPE 1


    //Parse the input file
    if(argc < 2) { printfln("Usage: %s inputfile_bosehubbardring",argv[0]); return 0; }
    auto input = InputGroup(argv[1],"input");
    string prefix="";
    if(argc > 2) {prefix=argv[2]; println("Use prefix: ",prefix);}

    auto N = input.getInt("N");
    auto Npart = input.getInt("Npart",N); //number of particles, default is N (half filling)
    auto Npartdown = input.getInt("Npartdown",N); //number of down particles, default is N (half filling)

    auto nsweeps = input.getInt("nsweeps");
    auto t1 = input.getReal("t1",1);
    auto t2 = input.getReal("t2",0);
    auto K = input.getReal("K",1); //lead coupling
    auto K2 = input.getReal("K2",1000); //Second lead coupling, set to 1000 to ignore and replace with value of K
    auto K3 = input.getReal("K3",0); //Third lead coupling for y-junction
    auto KL= input.getReal("KL",0); //coupling for ladder
    auto leadlength = input.getInt("leads",0); //lead length
    auto preleadlength = input.getInt("leadspre",0); //additional leadlength in front of first lead, Plead and Uleads apply only to that
    auto phi =input.getReal("phi",0); //has period 1
    auto phileads =input.getReal("phileads",0); //flux on leads has period 1
    auto U = input.getReal("U",0);
    auto Uleads = input.getReal("Uleads",0);
    auto V = input.getReal("V",0);
    auto Vcoupling=input.getReal("Vcoupling",1000);
    auto Vleads = input.getReal("Vleads",0);
    auto Pleads = input.getReal("Pleads",0);
    auto quiet = input.getYesNo("quiet",false);
    auto weightDMRG=input.getReal("weight",20.0);
    auto getdata=input.getInt("getdata",0);
    auto nstates = input.getInt("nstates",1);
    auto inittype=input.getInt("inittype",0);
    auto model=input.getInt("model",0);
    auto potential = input.getReal("potential",0);
    auto potential2 = input.getReal("potential2",0);
    auto dnpotentialmod = input.getReal("dnpotentialmod",1);
    auto sigmapotential = input.getReal("sigmapotential",2);
    auto centerpos=input.getInt("centerpos",0);


    auto phiEvolve =input.getInt("phiEvolve",0); // if 0, do not change phi during evolution, for 1, change at t=0 to phi=0


    auto impurity = input.getReal("impurity",0);
    auto impurity2 = input.getReal("impurity2",0);

    auto mpstype=input.getInt("mpstype",0);
    println("mpstype = ",mpstype);
    auto rangenp=input.getInt("rangenp",4);
    auto minnp=input.getInt("minnp",0);
    auto loadinitialoverlap=input.getInt("loadinitialoverlap",0);
    auto initialbondcompression=input.getInt("initialbondcompression",50);
    auto ttotal = input.getReal("ttotal",0.0);
    auto tau = input.getReal("tau",0.1);
    auto expectSnaps=input.getInt("expectSnaps",1);
    auto timecutoff = input.getReal("timecutoff",1E-13);
    auto timemaxm=input.getInt("timemaxm",500);
    auto timemethod=input.getInt("timemethod",0); //0: exact apply, 1: fit apply, 2: two step method exactapply, 3: two step method fitapply

    auto loadtime=input.getReal("loadtime",0); //Timestep from where to continue
    auto loadfile=input.getString("loadfile","0"); //Load filename if parameters changed
    auto stepsave=input.getInt("stepsave",0); //Save state at every step



	//Make dataset name to identify output
    auto dez=5;
    auto dataset="L"+to_string(N)+"N"+to_string(Npart)+"l"+to_string(leadlength)+"J"+ftos(t1,dez)+"g"+ftos(K,dez)+"f"+ftos(phi,dez)+"U"+ftos(U,dez)+"u"+ftos(Uleads,dez)+"V"+ftos(V,dez)+"v"+ftos(Vleads,dez)+"p"+ftos(Pleads,dez)+"P"+ftos(potential,dez);
    if(KL!=0) dataset+="K"+ftos(KL,dez);
    if(impurity!=0) dataset+=+"y"+ftos(impurity,dez);
    if(impurity2!=0) dataset+=+"Y"+ftos(impurity2,dez);
    if(K2!=1000) dataset+="k"+ftos(K2,dez); 
    if(K3!=0) dataset+="x"+ftos(K3,dez); 
    if(Vcoupling!=1000) dataset+="W"+ftos(Vcoupling,dez);
    dataset+="i"+to_string(inittype)+"d"+to_string(model)+"s"+to_string(nstates)+"T"+ftos(ttotal,dez)+"t"+ftos(tau,dez)+"m"+to_string(mpstype)+"M"+to_string(timemaxm)+"e"+to_string(timemethod)+"c"+ftos(-log10(timecutoff),1);
    if(rangenp!=4) dataset+="r"+to_string(rangenp);
    if(minnp!=0) dataset+="m"+to_string(minnp);
    if(phileads!=0) dataset+="l"+ftos(phileads,dez);
    if(potential2!=0) dataset+="p"+ftos(potential2,dez);
    if(sigmapotential!=2) dataset+="s"+ftos(sigmapotential,dez);
    if(centerpos!=0) dataset+="c"+to_string(centerpos);
    if(getdata!=0) dataset+="G"+to_string(getdata);
    if(phiEvolve!=0) dataset+="h"+to_string(phiEvolve); 
    if(preleadlength>0) dataset+="L"+to_string(preleadlength);

    auto datasetbare=dataset;
    datasetbare=ReplaceAll(datasetbare,".","_");
    if(stepsave==0) dataset+="A"+ftos(loadtime+ttotal,dez); 
    dataset=ReplaceAll(dataset,".","_");
 
    auto datasetSave=dataset;
    string datasetdelete;
    println("basedataset: ",dataset);


    double K2save=K2;
    double K3save=K3;
    double Vcouplingsave=Vcoupling;
    if(K2==1000)
	{
	K2=K; //Replace K2 by K if K2=1000
	println("Replace K2 with K=",K2);
	}
     else
	{
	println("K2=",K2);
	}
  
    if(K3==1000){ //If set to 1000 reuse value from K
	K3=K;
	println("Replaced K3 with K=",K3);
    }

    if(Vcoupling==1000){//If set to 1000 reuse value from Vleads
	Vcoupling=Vleads; //This may has to be refined if you use preleads and different V
    }
	println("Vcoupling set to ",Vcoupling);


    auto table = InputGroup(input,"sweeps");
    
    auto sweeps = Sweeps(nsweeps,table);
    std::complex<double> phase;
    std::complex<double> Kphase;
    std::complex<double> phaseLoc;
    std::complex<double> phasedmrgLoc;
    std::complex<double> imag=std::complex<double> (0.0,1.0);


    println(sweeps);

    //
    // Initialize the site degrees of freedom.
    //
    int Ntotal=0;
    if(model==0 or model==3) Ntotal=N+2*leadlength+preleadlength; //Ring
    if(model==1) Ntotal=2*N+leadlength+preleadlength; //Y
    if(model==2) Ntotal=N+leadlength+preleadlength; //Chain
    if(model==4) Ntotal=3*N+1+2*leadlength+preleadlength; //AB-Cage
    if(model==5) Ntotal=2*N+2*leadlength+preleadlength; //Ladder

    auto argsSiteset = Args();
    if(mpstype==7) argsSiteset=Args("rangenp=",rangenp,"minnp=",minnp);
    auto maxpsite=0;
    auto createOp="Adag"; //A stands for bosonic type operators. C will introduce fermionic anti-commuation rules on different sites when using AutoMPO by some internal rules of ITensor. For correlators, Jordan-Wigner has to be done by hand!  
    auto destroyOp="A";
    auto createOpdn="Cdagdn";
    auto destroyOpdn="Cdn";
    auto Nop="N";
    auto Nopdn="Ndn";
    auto NNop="NN";


    SiteSet sites;
    if(loadtime!=0)
    {
	string loadstring;
	if(loadfile=="0"){
		loadstring=datasetbare+"A"+ftos(loadtime,dez); //Construct file to load from
	}else{
		loadstring=loadfile; //Construct file to load from
	}
    	loadstring=ReplaceAll(loadstring,".","_");
	println("Reading site file from "+loadstring);
	loadsites(sites,"sites"+loadstring, mpstype,argsSiteset);
	/*if(mpstype==0){
		sites=readFromFile<BoseHubbard>("sites"+loadstring);
	}else if(mpstype==1 or mpstype==2){
		sites=readFromFile<Spinless>("sites"+loadstring);
	}else if(mpstype==3 or mpstype==6){
		sites=readFromFile<Hubbard>("sites"+loadstring);
	}else if(mpstype==4){
		sites=readFromFile<BoseHubbard3>("sites"+loadstring);
	}else if(mpstype==5){
		sites=readFromFile<BoseHubbard5>("sites"+loadstring);
	}*/
    }else
    {
	    if(mpstype==0) sites = BoseHubbard(Ntotal);
	    if(mpstype==1 or mpstype==2) sites = Spinless(Ntotal);
	    if(mpstype==3 or mpstype==6) sites = Hubbard(Ntotal);
	    if(mpstype==4) sites = BoseHubbard3(Ntotal);
	    if(mpstype==5) sites = BoseHubbard5(Ntotal);
	    if(mpstype==7) sites = BoseHubbardGen(Ntotal,argsSiteset); //Flexible atom number and background
    }



    if(mpstype==2)
    {
	createOp="Cdag";
	destroyOp="C";
    }

    if(mpstype==3 or mpstype==6)
    {
	if(mpstype==3){
		createOp="Cdagup";
		destroyOp="Cup";
	}
	if(mpstype==6){
		createOp="Adagup";
		destroyOp="Aup";
	    	createOpdn="Adagdn";
	    	destroyOpdn="Adn";
	}
	Nop="Nup";
	NNop="Nupdn";
    }


    if(mpstype==0) maxpsite=4;
    if(mpstype==4) maxpsite=3;
    if(mpstype==5) maxpsite=5;
    if(mpstype==1 or mpstype==2 or mpstype==3 or mpstype==6) maxpsite=1;
    if(mpstype==7) maxpsite=rangenp;

    println("Maxpsite:", maxpsite);

    //
    // Create the Hamiltonian using AutoMPO
    //
    //auto ampo = AutoMPO(sites);
    //auto ampoSource=AutoMPO(sites); //Special MPO for lead only initialization
    //auto ampoBump = AutoMPO(sites);



    auto Nmain=N;
    if(model==4) Nmain=3*N+1;
    if(model==5) Nmain=2*N;

    //Define variables for leads
    auto firstleadlength=leadlength+preleadlength; //first lead, is leadlength+preleadlength
    auto firststageleadlength=leadlength; //is either leadlength (without preleads), or preleadlength 
    if(preleadlength>0) firststageleadlength=preleadlength;//Apply Pleads, Uleads,Vleads only to the first preleadlength sites of source lead

    auto paramvalues=new double[1];
    auto energyList=new double[1];
    auto placeholderP=new double[1];

    double** dens = new double*[Ntotal];
    double** densdn = new double*[Ntotal];
    double** densdens = new double*[Ntotal];
    double** currentP = new double*[Ntotal];
    for(int i = 0; i < Ntotal; i++)
	{
	      dens[i] = new double[1];
	      densdn[i] = new double[1];
	      densdens[i] = new double[1];
	      currentP[i] = new double[1];
       }




 

   hblueprint hdmrg;
   hblueprint hsourceonly;

   double KLoc;
   double K2Loc;
   double t1Left;
   double t1Right;


   IQMPS psi;
   hblueprint hevolve;

   int Nsystem;

//Run DMRG to prepare ground state

    hdmrg =hblueprint();
    hsourceonly=hblueprint();

	//Set phases for Y-junction and rings in tunneling terms
    phase=exp(imag*2*M_PI*phi/N);

	
    Kphase=1;
    if(model==1 or model==2) phase=1;
    if(model==1 and K3!=0) {
	Kphase=exp(imag*2*M_PI*phi/3);
	println("Yjunction phase", Kphase);
    }
    if(model==4){
	 phase=exp(imag*2*M_PI*phi/4);
	println("AB cage phase",phase);
    }
    if(model==5){
	 phase=exp(imag*M_PI*phi);
	println("Ladder phase",phase);
    }

    println("complex phase coupling",phase);
    //For model==3, assign coupling in Y-junction fashion
    KLoc=K;
    if(model==3) KLoc=t1;
    K2Loc=K2;
    if(model==3) K2Loc=t1;
    t1Left=t1;
    if(model==3) t1Left=K;
    t1Right=t1;
    if(model==3) t1Right=K2;


	//Add all the entries for the Hamiltonian
    auto phaseleads=exp(imag*2*M_PI*phileads/firststageleadlength);//Phase applied to leads
    if(firstleadlength>0){
    	//Do lead L
	    for(int b = 1; b < firstleadlength; ++b)
		{
		if(b<firststageleadlength){
			hsourceonly.add(-t1*phaseleads,createOp,b,destroyOp,b+1,true);
		    	//ampoSource += -t1*phaseleads,createOp,b,destroyOp,b+1;
		    	//ampoSource += -t1*conj(phaseleads),createOp,b+1,destroyOp,b;
		}
		//Source only Hamiltonian only for first stage
		/*else{
			hsourceonly.add(-t1,createOp,b,destroyOp,b+1,true);
		    	//ampoSource += -t1,createOp,b,destroyOp,b+1;
		    	//ampoSource += -t1,createOp,b+1,destroyOp,b;
		}*/
		hdmrg.add(-t1*phaseleads,createOp,b,destroyOp,b+1,true);
	    	//ampo += -t1,createOp,b,destroyOp,b+1;
	    	//ampo += -t1,createOp,b+1,destroyOp,b;

		if(b<=firststageleadlength){
			hdmrg.add(Vleads,Nop,b,Nop,b+1,false);
			hsourceonly.add(Vleads,Nop,b,Nop,b+1,false);
			//if(Vleads!=0) ampo += Vleads,Nop,b,Nop,b+1;
			//if(Vleads!=0) ampoSource += Vleads,Nop,b,Nop,b+1;
		}else{
			hdmrg.add(V,Nop,b,Nop,b+1,false);
			hsourceonly.add(V,Nop,b,Nop,b+1,false);
			//if(V!=0) ampo += V,Nop,b,Nop,b+1;
			//if(V!=0) ampoSource += V,Nop,b,Nop,b+1;
		}
		if(mpstype==3 or mpstype==6)
			{
		    	//ampo += -t1,createOpdn,b,destroyOpdn,b+1;
		    	//ampo += -t1,createOpdn,b+1,destroyOpdn,b;
			hdmrg.add(-t1,createOpdn,b,destroyOpdn,b+1,true);
			if(b<firststageleadlength){
				hsourceonly.add(-t1*phaseleads,createOpdn,b,destroyOpdn,b+1,true);
			    	//ampoSource += -t1*phaseleads,createOpdn,b,destroyOpdn,b+1;
			    	//ampoSource += -t1*conj(phaseleads),createOpdn,b+1,destroyOpdn,b;
			}/*else{
				hsourceonly.add(-t1,createOpdn,b,destroyOpdn,b+1,true);
			    	//ampoSource += -t1,createOpdn,b,destroyOpdn,b+1;
			    	//ampoSource += -t1,createOpdn,b+1,destroyOpdn,b;
			}*/
			hdmrg.add(Vleads,Nopdn,b,Nopdn,b+1,false);
			hsourceonly.add(Vleads,Nopdn,b,Nopdn,b+1,false);
        		//if(Vleads!=0) ampo += Vleads,Nopdn,b,Nopdn,b+1;
        		//if(Vleads!=0) ampoSource += Vleads,Nopdn,b,Nopdn,b+1;
			}
		}


    	//Do ring-lead L coupling or first Y-junction coupling or first ladder coupling
	hdmrg.add(-KLoc*Kphase,createOp,firstleadlength,destroyOp,1+firstleadlength,true);
	hdmrg.add(Vcoupling,Nop,firstleadlength,Nop,1+firstleadlength,false);
    	//ampo += -KLoc*Kphase,createOp,firstleadlength,destroyOp,1+firstleadlength;
    	//ampo += -KLoc*conj(Kphase),createOp,1+firstleadlength,destroyOp,firstleadlength;
	/*
	if(preleadlength>0){
		if(V!=0) hdmrg.add(V,Nop,firstleadlength,Nop,1+firstleadlength);
        	//if(V!=0) ampo += V,Nop,firstleadlength,Nop,1+firstleadlength;
	}else{
		if(Vleads!=0) hdmrg.add(Vleads,Nop,firstleadlength,Nop,1+firstleadlength);
        	//if(Vleads!=0) ampo += Vleads,Nop,firstleadlength,Nop,1+firstleadlength;
	}
	*/
	if(mpstype==3 or mpstype==6)
	{
		hdmrg.add(-KLoc*Kphase,createOpdn,firstleadlength,destroyOpdn,1+firstleadlength,true);
		hdmrg.add( Vcoupling,Nopdn,firstleadlength,Nopdn,1+firstleadlength,false);
	    	//ampo += -KLoc*Kphase,createOpdn,firstleadlength,destroyOpdn,1+firstleadlength;
	    	//ampo += -KLoc*conj(Kphase),createOpdn,1+firstleadlength,destroyOpdn,firstleadlength;
		/*
		if(preleadlength>0){
			if(V!=0) hdmrg.add( V,Nopdn,firstleadlength,Nopdn,1+firstleadlength);
        		//if(V!=0) ampo += V,Nopdn,firstleadlength,Nopdn,1+firstleadlength;
		}else{
			if(Vleads!=0) hdmrg.add( Vleads,Nopdn,firstleadlength,Nopdn,1+firstleadlength);
        		//if(Vleads!=0) ampo += Vleads,Nopdn,firstleadlength,Nopdn,1+firstleadlength;
		}
		*/
	}

	if(model==5){
	    	//Do second coupling to ladder
		hdmrg.add(-KLoc*Kphase,createOp,firstleadlength,destroyOp,2+firstleadlength,true);
		hdmrg.add(Vcoupling,Nop,firstleadlength,Nop,2+firstleadlength,false);
		if(mpstype==3 or mpstype==6)
		{
			hdmrg.add(-KLoc*Kphase,createOpdn,firstleadlength,destroyOpdn,2+firstleadlength,true);
			hdmrg.add( Vcoupling,Nopdn,firstleadlength,Nopdn,2+firstleadlength,false);

		}

	}

	//Do second drain lead coupling for Y-junction model
	if(model==1)
	{
		hdmrg.add(-K2*conj(Kphase),createOp,firstleadlength,destroyOp,2+firstleadlength,true);
	    	//ampo += -K2*conj(Kphase),createOp,firstleadlength,destroyOp,2+firstleadlength;
	    	//ampo += -K2*Kphase,createOp,2+firstleadlength,destroyOp,firstleadlength;
		hdmrg.add(Vcoupling,Nop,firstleadlength,Nop,2+firstleadlength,false);
		if(mpstype==3 or mpstype==6)
		{
			hdmrg.add(-K2*conj(Kphase),createOpdn,firstleadlength,destroyOpdn,2+firstleadlength,true);
		    	//ampo += -K2*conj(Kphase),createOpdn,firstleadlength,destroyOpdn,2+firstleadlength;
		    	//ampo += -K2*Kphase,createOpdn,2+firstleadlength,destroyOpdn,firstleadlength;
			hdmrg.add(Vcoupling,Nopdn,firstleadlength,Nopdn,2+firstleadlength,false);
	
		}

	}

	//Do third drain lead coupling for Y-junction model
	if(model==1 and K3!=0)
	{
		hdmrg.add(-K3*Kphase,createOp,1+firstleadlength,destroyOp,2+firstleadlength,true);
	    	//ampo += -K3*Kphase,createOp,1+firstleadlength,destroyOp,2+firstleadlength;
	    	//ampo += -K3*conj(Kphase),createOp,2+firstleadlength,destroyOp,1+firstleadlength;
		hdmrg.add(Vcoupling,Nop,1+firstleadlength,Nop,2+firstleadlength,false);
		if(mpstype==3 or mpstype==6)
		{
			hdmrg.add(-K3*Kphase,createOpdn,1+firstleadlength,destroyOpdn,2+firstleadlength,true);
		    	//ampo += -K3*Kphase,createOpdn,1+firstleadlength,destroyOpdn,2+firstleadlength;
		    	//ampo += -K3*conj(Kphase),createOpdn,2+firstleadlength,destroyOpdn,1+firstleadlength;
			hdmrg.add(Vcoupling,Nopdn,1+firstleadlength,Nopdn,2+firstleadlength,false);
		}

	}

	if(Pleads!=0)
	{

	for(int i = 1; i <= firststageleadlength; ++i) 
		{
		hdmrg.add(Pleads,Nop,i);
		if ((model==0 or model==3 or model==4 or model==5) and preleadlength==0) hdmrg.add(Pleads,Nop,i+Nmain+firstleadlength);
		hsourceonly.add(Pleads,Nop,i);
		//ampo += Pleads,Nop,i;
		//if ((model==0 or model==3) and preleadlength==0) ampo += Pleads,Nop,i+N+firstleadlength;
		//ampoSource += Pleads,Nop,i;    //Do special MPO for leads only evolution
		if(mpstype==3 or mpstype==6)
		{
			hdmrg.add( Pleads,Nopdn,i);
			if ((model==0 or model==3 or model==4 or model==5) and preleadlength==0) hdmrg.add(Pleads,Nopdn,i+Nmain+firstleadlength);
			hsourceonly.add( Pleads,Nopdn,i);
			//ampo += Pleads,Nopdn,i;
			//if ((model==0 or model==3) and preleadlength==0) ampo += Pleads,Nopdn,i+N+firstleadlength;
			//ampoSource += Pleads,Nopdn,i;    //Do special MPO for leads only evolution
		}
		}

	}


    }


    if(maxpsite>1 or (mpstype==3 or mpstype==6))
	{
	    auto Nsystem=N;//Over what range to apply U
	    if (model==1) Nsystem=2*N;
	    if((model==1 or model==2) and preleadlength>0) Nsystem+=leadlength;
	    if((model==0 or model==3) and preleadlength>0) Nsystem=N+2*leadlength;
	    if (model==4) Nsystem=3*N+1;
	    if(model==4 and preleadlength>0) Nsystem=3*N+1+2*leadlength;
	    if (model==5) Nsystem=2*N;
	    if(model==5 and preleadlength>0) Nsystem=2*N+2*leadlength;
	    for(int i = 1; i <= Nsystem; ++i) 
		{
		if(U!=0)
		{
			if(mpstype==3 or mpstype==6)
				{
				hdmrg.add(U,"Nupdn",i+firststageleadlength);
				//ampo += U,"Nupdn",i+firststageleadlength;
				}
			else
				{
				hdmrg.add(U/2,"NInt",i+firststageleadlength);
				//ampo += U/2,"NInt",i+firststageleadlength;
				}
			}
		}

		//Lead interaction 

	    for(int i = 1; i <= firststageleadlength; ++i) 
		{
		if(Uleads!=0)
		{
			if(mpstype==3 or mpstype==6)
				{
				hdmrg.add(Uleads,"Nupdn",i);
				if ((model==0 or model==3 or model==4 or model==5) and preleadlength==0) hdmrg.add( Uleads,"Nupdn",i+Nmain+firstleadlength);
				hsourceonly.add(Uleads,"Nupdn",i);
				//ampo += Uleads,"Nupdn",i;
				//if ((model==0 or model==3) and preleadlength==0) ampo += Uleads,"Nupdn",i+N+firstleadlength;
				//ampoSource += Uleads,"Nupdn",i;    //Do special MPO for leads only evolution
				}
			else
				{
				hdmrg.add(Uleads/2,"NInt",i);
				if ((model==0 or model==3 or model==4 or model==5) and preleadlength==0) hdmrg.add(Uleads/2,"NInt",i+Nmain+firstleadlength);
				hsourceonly.add(Uleads/2,"NInt",i);
				//ampo += Uleads/2,"NInt",i;
				//if ((model==0 or model==3) and preleadlength==0) ampo += Uleads/2,"NInt",i+N+firstleadlength;
				//ampoSource += Uleads/2,"NInt",i;    //Do special MPO for leads only evolution
				}
		}
		}
       }


    


    //Impurity
    if((model==0 or model==3) and impurity!=0)
    {
	int impuritypos=firstleadlength+ceil(N/2);
	if((N/2)%2==1) impuritypos-=1;
	hdmrg.add(impurity,Nop,impuritypos);
	if(mpstype==3 or mpstype==6) hdmrg.add( impurity,Nopdn,impuritypos);
	//ampo += impurity,Nop,impuritypos;
	//if(mpstype==3 or mpstype==6) ampo += impurity,Nopdn,impuritypos;
	println("Placed impurity at ",impuritypos);
    }
    if((model==0 or model==3) and impurity2!=0)
    {
	int impuritypos=firstleadlength+ceil(N/2)+1;
	if((N/2)%2==1) impuritypos+=1;
	hdmrg.add(impurity2,Nop,impuritypos);
	if(mpstype==3 or mpstype==6) hdmrg.add( impurity2,Nopdn,impuritypos);
	//ampo += impurity2,Nop,impuritypos;
	//if(mpstype==3 or mpstype==6) ampo += impurity2,Nopdn,impuritypos;
	println("Placed impurity2 at ",impuritypos);
    }

    //Do ring-lead R coupling
    if(leadlength>0 and (model==0 or model==3 or model==4 or model==5)){
	    auto leadringsite=Nmain;
	    hdmrg.add(-K2Loc,createOp,firstleadlength+leadringsite,destroyOp,1+firstleadlength+Nmain,true);
 	    hdmrg.add(Vcoupling,Nop,firstleadlength+leadringsite,Nop,1+firstleadlength+Nmain,false);
	    //if(Vleads!=0) hdmrg.add(Vleads,Nop,firstleadlength+leadringsite,Nop,1+firstleadlength+Nmain);
	    //ampo += -K2Loc,createOp,firstleadlength+leadringsite,destroyOp,1+firstleadlength+N;
	    //ampo += -K2Loc,createOp,1+firstleadlength+N,destroyOp,firstleadlength+leadringsite;
	    //if(Vleads!=0) ampo += Vleads,Nop,firstleadlength+leadringsite,Nop,1+firstleadlength+N;
	    if(mpstype==3 or mpstype==6)
		  {
		    hdmrg.add( -K2Loc,createOpdn,firstleadlength+leadringsite,destroyOpdn,1+firstleadlength+Nmain,true);
	 	    hdmrg.add(Vcoupling,Nopdn,firstleadlength+leadringsite,Nopdn,1+firstleadlength+Nmain,false);
		    //ampo += -K2Loc,createOpdn,firstleadlength+leadringsite,destroyOpdn,1+firstleadlength+N;
		    //ampo += -K2Loc,createOpdn,1+firstleadlength+N,destroyOpdn,firstleadlength+leadringsite;
	    	    //if(Vleads!=0) ampo += Vleads,Nopdn,firstleadlength+leadringsite,Nopdn,1+firstleadlength+N;
		  }
	    if(model==5){//ladder second coupling to leads
		    hdmrg.add(-K2Loc,createOp,firstleadlength+leadringsite-1,destroyOp,1+firstleadlength+Nmain,true);
		    hdmrg.add(Vcoupling,Nop,firstleadlength+leadringsite-1,Nop,1+firstleadlength+Nmain,false);
		    if(mpstype==3 or mpstype==6)
			  {
			    hdmrg.add( -K2Loc,createOpdn,firstleadlength+leadringsite-1,destroyOpdn,1+firstleadlength+Nmain,true);
		 	    hdmrg.add(Vcoupling,Nopdn,firstleadlength+leadringsite-1,Nopdn,1+firstleadlength+Nmain,false);
			  }

	    }

	    //Do lead R
	    for(int b = 1; b < leadlength; ++b)
		{
		hdmrg.add(-t1,createOp,b+firstleadlength+Nmain,destroyOp,b+1+firstleadlength+Nmain,true);
		hdmrg.add(Vleads,Nop,b+firstleadlength+Nmain,Nop,b+1+firstleadlength+Nmain,false);
	    	//ampo += -t1,createOp,b+firstleadlength+N,destroyOp,b+1+firstleadlength+N;
	    	//ampo += -t1,createOp,b+1+firstleadlength+N,destroyOp,b+firstleadlength+N;
        	//if(Vleads!=0) ampo += Vleads,Nop,b+firstleadlength+N,Nop,b+1+firstleadlength+N;
		if(mpstype==3 or mpstype==6)
			{
			hdmrg.add(-t1,createOpdn,b+firstleadlength+Nmain,destroyOpdn,b+1+firstleadlength+Nmain,true);
			hdmrg.add(Vleads,Nopdn,b+firstleadlength+Nmain,Nopdn,b+1+firstleadlength+Nmain,false);
		    	//ampo += -t1,createOpdn,b+firstleadlength+N,destroyOpdn,b+1+firstleadlength+N;
		    	//ampo += -t1,createOpdn,b+1+firstleadlength+N,destroyOpdn,b+firstleadlength+N;
        		//if(Vleads!=0) ampo += Vleads,Nopdn,b+firstleadlength+N,Nopdn,b+1+firstleadlength+N;
			}

		}
    }


    hevolve=hblueprint(hdmrg); //create hevolve, up until now both hdmrg and hevolve were the same
    //Do system/phase related stuff
    auto phasedmrg=phase;
    if(phiEvolve==1) phase=1;
	//println("Create seperate ring part for initial Hamiltonian");
    	//ampoBump= AutoMPO(ampo);//Create MPO with a gaussian potential offset in the lead

 //To create ring in MPS, one needs to reorder it into nearest-neighbor hoppings.
    if(model==0 or model==3)//First ring hopping
    {
	    hdmrg.add(-t1Left*phasedmrg,createOp,1+firstleadlength,destroyOp,2+firstleadlength,true);
	    hdmrg.add(V,Nop,1+firstleadlength,Nop,2+firstleadlength,false);
	    hevolve.add(-t1Left*phase,createOp,1+firstleadlength,destroyOp,2+firstleadlength,true);
	    hevolve.add(V,Nop,1+firstleadlength,Nop,2+firstleadlength,false);
	    //ampoBump += -t1Left*phase,createOp,1+firstleadlength,destroyOp,2+firstleadlength;
	    //ampoBump += -t1Left*conj(phase),createOp,2+firstleadlength,destroyOp,1+firstleadlength;
	    //if(V!=0) ampoBump += V,Nop,1+firstleadlength,Nop,2+firstleadlength;
	   if(mpstype==3 or mpstype==6)
	   {
		hdmrg.add(-t1Left*phasedmrg,createOpdn,1+firstleadlength,destroyOpdn,2+firstleadlength,true);
		hdmrg.add(V,Nopdn,1+firstleadlength,Nopdn,2+firstleadlength,false);
		hevolve.add(-t1Left*phase,createOpdn,1+firstleadlength,destroyOpdn,2+firstleadlength,true);
		hevolve.add(V,Nopdn,1+firstleadlength,Nopdn,2+firstleadlength,false);
		//ampoBump += -t1Left*phase,createOpdn,1+firstleadlength,destroyOpdn,2+firstleadlength;
		//ampoBump += -t1Left*conj(phase),createOpdn,2+firstleadlength,destroyOpdn,1+firstleadlength;
		//if(V!=0) ampoBump += V,Nopdn,1+firstleadlength,Nopdn,2+firstleadlength;
	   }
    }

    Nsystem=N;
    if(model==1) Nsystem=2*N;
    if(model==2) Nsystem=N+1;
    if(model==4) Nsystem=3*N+1;
    if(model==5) Nsystem=2*N;
    
    auto internalhoplength=2; //Go by two to implement ring/junction/ladder
    if(model==2 or model==4) internalhoplength=1; //Regular chain 
    
    auto t1Array=new double[Nsystem-2];
    for(int b = 1; b < Nsystem-1; ++b) 
    {
	t1Array[b-1]=t1;
    }
    if(model==3)
    {
	t1Array[0]=t1Left;
	t1Array[Nsystem-3]=t1Right;
    }


    if(model==0 or model==1 or model==2 or model==3 or model==5){
	    for(int b = 1; b < Nsystem-1; ++b)
		{
		if(b%2==1){
			phaseLoc=conj(phase);
			phasedmrgLoc=conj(phasedmrg);
		}else{
			phaseLoc=phase;
			phasedmrgLoc=phasedmrg;
		}
		if(model==5){//For ladder conjugated is correct phase
			phaseLoc=conj(phaseLoc);
			phasedmrgLoc=conj(phasedmrgLoc);
		}
		hdmrg.add(-t1Array[b-1]*phasedmrgLoc,createOp,b+firstleadlength,destroyOp,b+internalhoplength+firstleadlength,true);
		hdmrg.add(V,Nop,b+firstleadlength,Nop,b+internalhoplength+firstleadlength,false);
		hevolve.add(-t1Array[b-1]*phaseLoc,createOp,b+firstleadlength,destroyOp,b+internalhoplength+firstleadlength,true);
		hevolve.add(V,Nop,b+firstleadlength,Nop,b+internalhoplength+firstleadlength,false);
		//ampoBump += -t1Array[b-1]*phaseLoc,createOp,b+firstleadlength,destroyOp,b+internalhoplength+firstleadlength;
		//ampoBump += -t1Array[b-1]*conj(phaseLoc),createOp,b+internalhoplength+firstleadlength,destroyOp,b+firstleadlength;
		//if(V!=0) ampoBump += V,Nop,b+firstleadlength,Nop,b+internalhoplength+firstleadlength;
		if(mpstype==3 or mpstype==6)
		  {
			hdmrg.add(-t1Array[b-1]*phasedmrgLoc,createOpdn,b+firstleadlength,destroyOpdn,b+internalhoplength+firstleadlength,true);
			hdmrg.add(V,Nopdn,b+firstleadlength,Nopdn,b+internalhoplength+firstleadlength,false);
			hevolve.add(-t1Array[b-1]*phaseLoc,createOpdn,b+firstleadlength,destroyOpdn,b+internalhoplength+firstleadlength,true);
			hevolve.add(V,Nopdn,b+firstleadlength,Nopdn,b+internalhoplength+firstleadlength,false);
			//ampoBump += -t1Array[b-1]*phaseLoc,createOpdn,b+firstleadlength,destroyOpdn,b+internalhoplength+firstleadlength;
			//ampoBump += -t1Array[b-1]*conj(phaseLoc),createOpdn,b+internalhoplength+firstleadlength,destroyOpdn,b+firstleadlength;
			//if(V!=0) ampoBump += V,Nopdn,b+firstleadlength,Nopdn,b+internalhoplength+firstleadlength;
		  }
		}
    }

    if(model==4)
    {
	for(int b = 0; b < N; ++b)
	{
		for(int c=0; c<2;++c){
			if(c%2==1){
				phaseLoc=conj(phase);
				phasedmrgLoc=conj(phasedmrg);
			}else{
				phaseLoc=phase;
				phasedmrgLoc=phasedmrg;
			}
			hdmrg.add(-t1*phasedmrgLoc,createOp,3*b+1+2*c+firstleadlength,destroyOp,3*b+2+2*c+firstleadlength,true);
			hdmrg.add(V,Nop,3*b+1+2*c+firstleadlength,Nop,3*b+2+2*c+firstleadlength,false);
			hevolve.add(-t1*phaseLoc,createOp,3*b+1+2*c+firstleadlength,destroyOp,3*b+2+2*c+firstleadlength,true);
			hevolve.add(V,Nop,3*b+1+2*c+firstleadlength,Nop,3*b+2+2*c+firstleadlength,false);
			if(mpstype==3 or mpstype==6)
			  {
				hdmrg.add(-t1*phasedmrgLoc,createOpdn,3*b+1+2*c+firstleadlength,destroyOpdn,3*b+2+2*c+firstleadlength,true);
				hdmrg.add(V,Nop,3*b+1+2*c+firstleadlength,Nopdn,3*b+2+2*c+firstleadlength,false);
				hevolve.add(-t1*phaseLoc,createOpdn,3*b+1+2*c+firstleadlength,destroyOpdn,3*b+2+2*c+firstleadlength,true);
				hevolve.add(V,Nop,3*b+1+2*c+firstleadlength,Nopdn,3*b+2+2*c+firstleadlength,false);
			}

		}
		for(int c=0; c<2;++c){
			if(c%2==0){
				phaseLoc=conj(phase);
				phasedmrgLoc=conj(phasedmrg);
			}else{
				phaseLoc=phase;
				phasedmrgLoc=phasedmrg;
			}
			hdmrg.add(-t1*phasedmrgLoc,createOp,3*b+1+c+firstleadlength,destroyOp,3*b+1+c+2+firstleadlength,true);
			hdmrg.add(V,Nop,3*b+1+c+firstleadlength,Nop,3*b+1+c+2+firstleadlength,false);
			hevolve.add(-t1*phaseLoc,createOp,3*b+1+c+firstleadlength,destroyOp,3*b+1+c+2+firstleadlength,true);
			hevolve.add(V,Nop,3*b+1+c+firstleadlength,Nop,3*b+1+c+2+firstleadlength,false);
			if(mpstype==3 or mpstype==6)
			  {
				hdmrg.add(-t1*phasedmrgLoc,createOpdn,3*b+1+c+firstleadlength,destroyOpdn,3*b+1+c+2+firstleadlength,true);
				hdmrg.add(V,Nop,3*b+1+c+firstleadlength,Nopdn,3*b+1+c+2+firstleadlength,false);
				hevolve.add(-t1*phaseLoc,createOpdn,3*b+1+c+firstleadlength,destroyOpdn,3*b+1+c+2+firstleadlength,true);
				hevolve.add(V,Nop,3*b+1+c+firstleadlength,Nopdn,3*b+1+c+2+firstleadlength,false);

			}
		}

	}

    }
    //Ladder KL coupling between ladders
    if(model==5){
	    for(int b = 1; b <= N; ++b)
		{
		hdmrg.add(-KL,createOp,2*b-1+firstleadlength,destroyOp,2*b+firstleadlength,true);
		hdmrg.add(V,Nop,2*b-1+firstleadlength,Nop,2*b+firstleadlength,false);
		hevolve.add(-KL,createOp,2*b-1+firstleadlength,destroyOp,2*b+firstleadlength,true);
		hevolve.add(V,Nop,2*b-1+firstleadlength,Nop,2*b+firstleadlength,false);
		if(mpstype==3 or mpstype==6)
		  {
			hdmrg.add(-KL,createOpdn,2*b-1+firstleadlength,destroyOpdn,2*b+firstleadlength,true);
			hdmrg.add(V,Nopdn,2*b-1+firstleadlength,Nopdn,2*b+internalhoplength+firstleadlength,false);
			hevolve.add(-KL,createOpdn,2*b-1+firstleadlength,destroyOpdn,2*b+firstleadlength,true);
			hevolve.add(V,Nopdn,2*b-1+firstleadlength,Nopdn,2*b+internalhoplength+firstleadlength,false);
		  }
		}
    }


    //Close ring hopping
    if(model==0 or model==3)
    {
	    if(N%2==0){
		phaseLoc=conj(phase);
		phasedmrgLoc=conj(phasedmrg);
	    }else{
	    	phaseLoc=phase;
	    	phasedmrgLoc=phasedmrg;
	    }
	    hdmrg.add(-t1Right*phasedmrgLoc,createOp,N-1+firstleadlength,destroyOp,N+firstleadlength,true);
	    hdmrg.add(V,Nop,N-1+firstleadlength,Nop,N+firstleadlength,false);
	    hevolve.add(-t1Right*phaseLoc,createOp,N-1+firstleadlength,destroyOp,N+firstleadlength,true);
	    hevolve.add(V,Nop,N-1+firstleadlength,Nop,N+firstleadlength,false);
	    //ampoBump += -t1Right*phaseLoc,createOp,N-1+firstleadlength,destroyOp,N+firstleadlength;
	    //ampoBump += -t1Right*conj(phaseLoc),createOp,N+firstleadlength,destroyOp,N-1+firstleadlength;
	    //if(V!=0) ampoBump += V,Nop,N-1+firstleadlength,Nop,N+firstleadlength;
	    if(mpstype==3 or mpstype==6)
		  {
		    hdmrg.add(-t1Right*phasedmrgLoc,createOpdn,N-1+firstleadlength,destroyOpdn,N+firstleadlength,true);
		    hdmrg.add(V,Nopdn,N-1+firstleadlength,Nopdn,N+firstleadlength,false);
		    hevolve.add(-t1Right*phaseLoc,createOpdn,N-1+firstleadlength,destroyOpdn,N+firstleadlength,true);
		    hevolve.add(V,Nopdn,N-1+firstleadlength,Nopdn,N+firstleadlength,false);
		    //ampoBump += -t1Right*phaseLoc,createOpdn,N-1+firstleadlength,destroyOpdn,N+firstleadlength;
		    //ampoBump += -t1Right*conj(phaseLoc),createOpdn,N+firstleadlength,destroyOpdn,N-1+firstleadlength;
		    //if(V!=0) ampoBump += V,Nopdn,N-1+firstleadlength,Nopdn,N+firstleadlength;
		}
   }
 
    

	//order sites for ring in order to turn ring into next-nearest coupling chain
    //Only for ring models used
    auto correctOrderSites=new int[N]; // correctOrderSites[i] gives the correct physical ordering of MPS site i (start counting at 0)
    correctOrderSites[0]=0;
    correctOrderSites[N-1]=int(N/2);
    for(int i = 1; i < int(N/2); ++i) 
    {
	correctOrderSites[2*i-1]=i;
	correctOrderSites[2*i]=N-i;
    }

    

	//Add potential for bump, only present in the dmrg Hamiltonian, but no in the evolving one
    if(potential!=0)
    {
	    if(model==0 and firstleadlength==0) 
		{
		    double gausscenter=int(floor(N/2.0))+1+centerpos;
		    for(int i = 1; i <= N; ++i) 
		    {
			hdmrg.add(-potential*exp(-pow(1.0*(correctOrderSites[i-1]+1)-gausscenter,2.0)/(2*pow(sigmapotential,2.0))),Nop,i);
		    	//ampoBump += -potential*exp(-pow(1.0*(correctOrderSites[i-1]+1)-gausscenter,2.0)/(2*pow(sigmapotential,2.0))),Nop,i;
			if(mpstype==3 or mpstype==6)
				{
				hdmrg.add(-potential*exp(-pow(1.0*(correctOrderSites[i-1]+1)-gausscenter,2.0)/(2*pow(sigmapotential,2.0))),Nopdn,i);
			    	//ampoBump += -potential*dnpotentialmod*exp(-pow(1.0*(correctOrderSites[i-1]+1)-gausscenter,2.0)/(2*pow(sigmapotential,2.0))),Nopdn,i;
				}
		    }
			if(potential2!=0){
			    double gausscenter2=int(1+centerpos);
			    for(int i = 1; i <= N; ++i) 
			    {
				double dist=min(1.0*(correctOrderSites[i-1]+1)-gausscenter2,N-abs(1.0*(correctOrderSites[i-1]+1)-gausscenter2));
				hdmrg.add(-potential2*exp(-pow(dist,2.0)/(2*pow(sigmapotential,2.0))),Nop,i);
			    	//ampoBump += -potential2*exp(-pow(dist,2.0)/(2*pow(sigmapotential,2.0))),Nop,i;
				if(mpstype==3 or mpstype==6)
					{
					hdmrg.add(-potential2*dnpotentialmod*exp(-pow(dist,2.0)/(2*pow(sigmapotential,2.0))),Nopdn,i);
				    	//ampoBump += -potential2*dnpotentialmod*exp(-pow(dist,2.0)/(2*pow(sigmapotential,2.0))),Nopdn,i;
					}
			    }
			}
		}else
		{
		    double gausscenter=int(floor(firststageleadlength/2.0))+1+centerpos;
		    for(int i = 1; i <= firststageleadlength; ++i) 
		    {
			hdmrg.add(-potential*exp(-pow(1.0*i-gausscenter,2.0)/(2*pow(sigmapotential,2.0))),Nop,i);
			hsourceonly.add(-potential*exp(-pow(1.0*i-gausscenter,2.0)/(2*pow(sigmapotential,2.0))),Nop,i);
		    	//ampoBump += -potential*exp(-pow(1.0*i-gausscenter,2.0)/(2*pow(sigmapotential,2.0))),Nop,i;
		    	//ampoSource += -potential*exp(-pow(1.0*i-gausscenter,2.0)/(2*pow(sigmapotential,2.0))),Nop,i;
			if(mpstype==3 or mpstype==6)
				{
				hdmrg.add( -potential*dnpotentialmod*exp(-pow(1.0*i-gausscenter,2.0)/(2*pow(sigmapotential,2.0))),Nopdn,i);
				hsourceonly.add( -potential*dnpotentialmod*exp(-pow(1.0*i-gausscenter,2.0)/(2*pow(sigmapotential,2.0))),Nopdn,i);
			    	//ampoBump += -potential*dnpotentialmod*exp(-pow(1.0*i-gausscenter,2.0)/(2*pow(sigmapotential,2.0))),Nopdn,i;
			    	//ampoSource += -potential*dnpotentialmod*exp(-pow(1.0*i-gausscenter,2.0)/(2*pow(sigmapotential,2.0))),Nopdn,i;
				}

		    }
		}
    }


    //auto H = IQMPO(ampo);
    IQMPO H=hevolve.generateMPO(sites);

    

	if(inittype==1){
	println("DMRG source only");
	hsourceonly.print();
	}else{
		println("DMRG");
		hdmrg.print();
	}
	println("Evolve");
	hevolve.print();
    
    if(loadtime!=0)
    {
	string loadstring;
	if(loadfile=="0"){
		loadstring=datasetbare+"A"+ftos(loadtime,dez); //Construct file to load from
	}else{
		loadstring=loadfile; //Construct file to load from
	}
    	loadstring=ReplaceAll(loadstring,".","_");
	println("Starting reading from "+loadstring);
	//readFromFile("sites"+loadstring,sites);
	//state = InitState(sitesA);
	//psi = IQMPS(state);
	psi=readFromFile<IQMPS>("psi"+loadstring,sites);
        println("Loaded state from "+loadstring);
    	Print(totalQN(psi));

    }
    else//begin load else
    {

    //
    // Set the initial wavefunction matrix product state
	//distribute particles across the system, in a random fashion that minimses the effort for dmrg to find ground state
    //
    auto state = InitState(sites);
    int p = Npart;
    int pdown=Npartdown;
    auto sitepnumber= new int[Ntotal];
    for(int i =0; i<  Ntotal; ++i) 
	{
	sitepnumber[i]=0;
	}

    auto NsitesDistribute=Ntotal;
    if(inittype==1) NsitesDistribute=firststageleadlength;
    auto setparticle = int(ceil(1.0*Npart/NsitesDistribute));
    if(setparticle>maxpsite) throw std::invalid_argument( "Number of particles too large for lattice" );
    for(int j=1;j<=setparticle+1;++j)
	{
    	    for(int i =0; i<  NsitesDistribute; ++i) 
		{
		    if(p>0)
			{
		    	sitepnumber[i]+=1;
			p-=1;
			}
		}
	}

    random_shuffle(&sitepnumber[0], &sitepnumber[NsitesDistribute-1]);
    println("Sites distribution:");
    for(int i =0; i<  Ntotal; ++i)
	{
	cout << sitepnumber[i]<< " ";
        }
    println("");
    for(int i =0; i<  Ntotal; ++i)
	{
	    auto s = std::to_string(sitepnumber[i]);
	    if((mpstype==3 or mpstype==6) and sitepnumber[i]>0) s="S"; //Currently assuming equal occupation for spin-fermions
	    state.set(i+1,s);
	}





    psi = IQMPS(state);

    Print(totalQN(psi));
    
    //
    // Begin the DMRG calculation
    //

    auto energy=0.0;
    IQMPO HBump;
    if(nstates>0)
     {
	println("Start DMRG");
	if(inittype==1)
		{
    		//auto Hleads = IQMPO(ampoSource);
		auto Hleads=hsourceonly.generateMPO(sites);
	    	energy = dmrg(psi,Hleads,sweeps,{"Quiet",quiet});
		}
	else
	if(inittype==2)
		{
    		//HBump = IQMPO(ampoBump);
		HBump=hdmrg.generateMPO(sites);
		println("Use HBump");
	    	energy = dmrg(psi,HBump,sweeps,{"Quiet",quiet});
		}
	else
		{
		energy = dmrg(psi,H,sweeps,{"Quiet",quiet});
		}

	if(initialbondcompression>0){
	    	auto args = Args("Cutoff=",timecutoff,"Maxm=",initialbondcompression);
		psi.orthogonalize(args);

		auto normcompmress=normalize(psi);
		printfln("Bond size, norm after compression= %d, %.5f",maxM(psi),normcompmress);
	}
    }

	//Calculate energy after dmrg
    if(inittype==2)
    {
    	energy = overlap(psi,HBump,psi);
    }
    else{
    	energy = overlap(psi,H,psi);
    }

    energyList[0]=energy;
    placeholderP[0]=0;
    //
    // Print the final energy reported by DMRG
    //
    printfln("\nGround State Energy = %.10f",energy);

    //
    // Measure densities
    //

    /*
    auto center = (int)(N/2);
    psi.position(center); 

    IQTensor wf = psi.A(center)*psi.A(center+1);
    auto Um = psi.A(center);
    IQTensor Sm,Vm;
    auto spectrum = svd(wf,Um,Sm,Vm);

    auto SvN = 0.;
    for(auto k : spectrum.eigs())
        {
        if(k > 1E-12) SvN += -k*log(k);
        }
    printfln("Across bond b=%d, SvN = %.10f",center,SvN);
    */



	//Copy and pasted from evolution section, may not respect quench related stuff for phase...
	auto poscurrenta=firstleadlength+1;
	auto poscurrentb=firstleadlength+2;
	auto currentcoupling=t1Left;
	auto phase1=phase;
	auto phase2=conj(phase);
	auto poscurrent2a=firstleadlength+1;
	auto poscurrent2b=firstleadlength+3;
	auto currentcoupling2=t1Left;
	auto phase3=conj(phase);
	auto poscurrent3a=firstleadlength+1;
	auto poscurrent3b=firstleadlength+3;
	auto currentcoupling3=t1Left;
	if(model==2 or model==1){
		phase1=Kphase;
	 	currentcoupling=K;
		poscurrenta=firstleadlength;
		poscurrentb=firstleadlength+1;
		phase2=conj(Kphase);
		poscurrent2a=firstleadlength;
		poscurrent2b=firstleadlength+2;
		currentcoupling2=K2;
		phase3=Kphase;
		poscurrent3a=firstleadlength+1;
		poscurrent3b=firstleadlength+2;
		currentcoupling3=K3;
	}
	if(model==4){
		phase1=1;
	 	currentcoupling=K;
		poscurrenta=firstleadlength;
		poscurrentb=firstleadlength+1;
		phase2=phase;
		poscurrent2a=firstleadlength+3*int(N/2)+1;
		poscurrent2b=firstleadlength+3*int(N/2)+2;
		currentcoupling2=t1;
		phase3=phase;
		poscurrent3a=firstleadlength+3*int(N/2)+2;
		poscurrent3b=firstleadlength+3*int(N/2)+4;
		currentcoupling3=t1;
	}
	if(model==5){
		phase1=phase;
		poscurrenta=firstleadlength+2*int(N/2)+1;
		poscurrentb=firstleadlength+2*int(N/2)+3;
		currentcoupling=t1;
		phase2=conj(phase);
	 	currentcoupling2=t1;
		poscurrent2a=firstleadlength+2*int(N/2)+2;
		poscurrent2b=firstleadlength+2*int(N/2)+4;
		phase3=1;
		poscurrent3a=firstleadlength+2*int(N/2)+1;
		poscurrent3b=firstleadlength+2*int(N/2)+2;
		currentcoupling3=KL;
	}

	
	auto param_index=0; //Unused index, just keep at zero
		//Calculate currents of dmrg state
	currentP[0][param_index]=-2*(-currentcoupling*phase1*getcorrelatorOverlap(poscurrenta,poscurrentb,createOp,destroyOp,sites,psi,psi)).imag(); //Calculates Adag(a)*a(b), 
	currentP[1][param_index]=-2*(-currentcoupling2*phase2*getcorrelatorOverlap(poscurrent2a,poscurrent2b,createOp,destroyOp,sites,psi,psi)).imag(); //Calculates Adag(a)*a(b), 
	currentP[2][param_index]=-2*(-currentcoupling3*phase3*getcorrelatorOverlap(poscurrent3a,poscurrent3b,createOp,destroyOp,sites,psi,psi)).imag(); //Calculates Adag(a)*a(b), 

   
    //Vector dens(Ntotal),densdn(Ntotal),densdens(Ntotal);
    for(int j = 1; j <= Ntotal; ++j)
        {
        psi.position(j);
        //dens(j-1) = (dag(prime(psi.A(j),Site))*sites.op(Nop,j)*psi.A(j)).real();
	dens[j-1][param_index] = getlocalExpectation(j, Nop,sites, psi);
	if(maxpsite>1 or (mpstype==3 or mpstype==6))
		{
			//densdens(j-1) = (dag(prime(psi.A(j),Site))*sites.op(NNop,j)*psi.A(j)).real();
			densdens[j-1][param_index] = getlocalExpectation(j, NNop,sites, psi);
			if(mpstype==3 or mpstype==6)
			{
				densdn[j-1][param_index] = getlocalExpectation(j, Nopdn,sites, psi);
				//densdn(j-1)=(dag(prime(psi.A(j),Site))*sites.op(Nopdn,j)*psi.A(j)).real();
			}
		}
		else
		{
			densdens[j-1][param_index] = dens[j-1][param_index];
		}
        }

    
	    println("Density:");
	    for(int j = 0; j < Ntotal; ++j)
		{
		if(mpstype==3 or mpstype==6)
			{
			printfln("%d %.6f, %.6f, %.6f",1+j,dens[j][param_index],densdn[j][param_index],densdens[j][param_index]-dens[j][param_index]*densdn[j][param_index]);
			}
		else
			{
			printfln("%d %.6f, %.6f",1+j,dens[j][param_index],densdens[j][param_index]-dens[j][param_index]*dens[j][param_index]);
			}
		}


	    auto totalpn=0.0;
	    auto totalsourcepn=0.0;
	    auto totaldrainpn=0.0;
	    auto totalpreleadpn=0.0;
	    auto totalpostleadpn=0.0;
	    for(int j=0;j<Ntotal;++j)
		totalpn+=dens[j][param_index];
	    for(int j=0;j<firstleadlength;++j)
		totalsourcepn+=dens[j][param_index];
	    if(preleadlength>0) {
		    for(int j=0;j<firststageleadlength;++j)
			totalpreleadpn+=dens[j][param_index];
		    for(int j=0;j<leadlength;++j)
			totalpostleadpn+=dens[j+firststageleadlength][param_index];
	    }
	    if(model==1){
		    for(int j=0;j<2*N;++j){
			totaldrainpn+=dens[-j-1+Ntotal][param_index];
			}
	    }else if(model==2){
		    for(int j=0;j<N;++j){
			totaldrainpn+=dens[-j-1+Ntotal][param_index];
			}
	    }else if(model==0 or model==3 or model==4 or model==5){
		    for(int j=0;j<leadlength;++j){
			totaldrainpn+=dens[-j-1+Ntotal][param_index];
		    }

	    }

	    println("Total particle number ", totalpn);
	    println("Total particle number source: ", totalsourcepn,", Total particle number drain: ", totaldrainpn);
	    if(preleadlength>0) println("Pnumber prelead: ", totalpreleadpn,", Pnumber postlead: ", totalpostleadpn);

	    if(mpstype==3 or mpstype==6)
		{
		auto diffupdn=0.0;
		for(int j=0;j<Ntotal;++j){
			diffupdn+=abs(dens[j][param_index]-densdn[j][param_index]);
		}
		println("Difference up/dn ", diffupdn);
		}
		
	}//End load else

   

	//DMRG related stuff finished


	//Do time evolution
    double finishtime=0.0;
    //Time evolution
    if(ttotal>0)
	{
	IQMPS initpsi;
	if(loadinitialoverlap==1){
		string datasetInit;
		if(loadtime==0){
			//Write initial wavefunction to file

			datasetInit=dataset+"A"+ftos(0,dez);
			datasetInit=ReplaceAll(datasetInit,".","_");
			println("Write init wavefunction to file",datasetInit);
			writeToFile("psi"+datasetInit,psi);
			writeToFile("sites"+datasetInit,sites); //Write site object
			initpsi=IQMPS(psi);

		}else{

			SiteSet initsites;
			datasetInit=dataset+"A"+ftos(0,dez);
			datasetInit=ReplaceAll(datasetInit,".","_");
			println("Read initial wavefunction from file", datasetInit);
			loadsites(initsites,"sites"+datasetInit, mpstype,argsSiteset);
			initpsi=readFromFile<IQMPS>("psi"+datasetInit,sites);
			println("Succesfully loaded sites and wavefunction");

		}
	}else{
		initpsi=IQMPS(psi);
	}


	std::complex<double> taua;
	std::complex<double> taub;
	if(timemethod==0 or timemethod==1) 
	{
		taua=tau*Cplx_i;
		taub=0;

	}else
	if(timemethod==2 or timemethod==3)  //Combine two timesteps for tau^3 accuracy as described in arXiv 1407.1832
	{
		taua=(1+Cplx_i)/2*Cplx_i*tau;
		taub=(1-Cplx_i)/2*Cplx_i*tau;
	}

	//auto expH = toExpH<IQTensor>(ampo,taua);
	//auto expHb = toExpH<IQTensor>(ampo,taub);
	auto expH= hevolve.generateMPO(sites,taua);
	auto expHb= hevolve.generateMPO(sites,taub);

	auto args = Args("Cutoff=",timecutoff,"Maxm=",timemaxm,"Verbose=",false);
	auto nt = int(ttotal/tau+(1e-9*(ttotal/tau)));

	//exactApplyMPO(psi,expH,psi);

	
	auto numberSnaps=int(nt/expectSnaps);
	auto timearray= new double[numberSnaps+1];
	//auto entropy= new double[numberSnaps+1];
	auto nrm= new double[numberSnaps+1];
	auto bond= new int[numberSnaps+1];


	auto comptime= new double[numberSnaps+1];
	auto placeholder=new double[numberSnaps+1];
	auto overlapinit=new double[numberSnaps+1];
	timearray[0]=loadtime;

	for(int n = 1; n <=numberSnaps; ++n)
	{
		timearray[n]=loadtime+n*tau*expectSnaps;
	}
	finishtime=timearray[numberSnaps];
	double** densR = new double*[Ntotal];
	double** densdnR = new double*[Ntotal];
	double** densdensR = new double*[Ntotal];
	double** current = new double*[Ntotal];

	  for(int i = 0; i < Ntotal; i++)
		{
	      densR[i] = new double[numberSnaps+1];
	      densdnR[i] = new double[numberSnaps+1];
	      densdensR[i] = new double[numberSnaps+1];
	      current[i] = new double[numberSnaps+1];

		}
	nrm[0]=normalize(psi);
	
	//Set all the parameters again, here for time evolution#
	//Ring-lead
	auto poscurrenta=firstleadlength+1;
	auto poscurrentb=firstleadlength+2;
	auto currentcoupling=t1Left;
	auto phase1=phase;
	auto phase2=conj(phase);
	auto poscurrent2a=firstleadlength+1;
	auto poscurrent2b=firstleadlength+3;
	auto currentcoupling2=t1Left;
	auto phase3=conj(phase);
	auto poscurrent3a=firstleadlength+1;
	auto poscurrent3b=firstleadlength+3;
	auto currentcoupling3=t1Left;
	
	
	if(model==2 or model==1){ //Y-junction
		phase1=Kphase;
	 	currentcoupling=K;
		poscurrenta=firstleadlength;
		poscurrentb=firstleadlength+1;
		phase2=conj(Kphase);
		poscurrent2a=firstleadlength;
		poscurrent2b=firstleadlength+2;
		currentcoupling2=K2;
		phase3=Kphase;
		poscurrent3a=firstleadlength+1;
		poscurrent3b=firstleadlength+2;
		currentcoupling3=K3;
	}
	if(model==4){
		phase1=1;
	 	currentcoupling=K;
		poscurrenta=firstleadlength;
		poscurrentb=firstleadlength+1;
		phase2=phase;
		poscurrent2a=firstleadlength+3*int(N/2)+1;
		poscurrent2b=firstleadlength+3*int(N/2)+2;
		currentcoupling2=t1;
		phase3=phase;
		poscurrent3a=firstleadlength+3*int(N/2)+2;
		poscurrent3b=firstleadlength+3*int(N/2)+4;
		currentcoupling3=t1;
	}
	if(model==5){
		phase1=phase;
		poscurrenta=firstleadlength+2*int(N/2)+1;
		poscurrentb=firstleadlength+2*int(N/2)+3;
		currentcoupling=t1;
		phase2=conj(phase);
	 	currentcoupling2=t1;
		poscurrent2a=firstleadlength+2*int(N/2)+2;
		poscurrent2b=firstleadlength+2*int(N/2)+4;
		phase3=1;
		poscurrent3a=firstleadlength+2*int(N/2)+1;
		poscurrent3b=firstleadlength+2*int(N/2)+2;
		currentcoupling3=KL;
	}



	current[0][0]=-2*(-currentcoupling*phase1*getcorrelatorOverlap(poscurrenta,poscurrentb,createOp,destroyOp,sites,psi,psi)).imag(); //Calculates Adag(a)*a(b), 
	current[1][0]=-2*(-currentcoupling2*phase2*getcorrelatorOverlap(poscurrent2a,poscurrent2b,createOp,destroyOp,sites,psi,psi)).imag(); //Calculates Adag(a)*a(b), 
	current[2][0]=-2*(-currentcoupling3*phase3*getcorrelatorOverlap(poscurrent3a,poscurrent3b,createOp,destroyOp,sites,psi,psi)).imag(); //Calculates Adag(a)*a(b), 
		println("Current t=0 :",current[0][0],", ",current[1][0],", ",current[2][0]);

	overlapinit[0]=pow(abs(overlapC(initpsi,psi)),2);
	for(int j = 1; j <= Ntotal; ++j)
	{
		psi.position(j);
		densR[j-1][0] = getlocalExpectation(j, Nop,sites, psi);
		//densR[j-1][0] = (dag(prime(psi.A(j),Site))*sites.op(Nop,j)*psi.A(j)).real();
		if(mpstype==3 or mpstype==6)
			{
			//densdnR[j-1][0] = (dag(prime(psi.A(j),Site))*sites.op(Nopdn,j)*psi.A(j)).real();
			//densdensR[j-1][0] = (dag(prime(psi.A(j),Site))*sites.op(NNop,j)*psi.A(j)).real();
			densdnR[j-1][0] = getlocalExpectation(j, Nopdn,sites, psi);

			}
		if(maxpsite>1 or (mpstype==3 or mpstype==6)){
			densdensR[j-1][0] = getlocalExpectation(j, NNop,sites, psi);
		}

	}


	bond[0]=maxM(psi);
    	comptime[0] = ( std::clock() - start ) / (double) CLOCKS_PER_SEC; //Time until start of time evolution
	placeholder[0]=0;
        printfln("\nComputation time until now= %.4f",comptime[0]);
	for(int n = 1; n <= numberSnaps; ++n)
	    {
		auto currentnumberSnaps=n;
		//Normalize psi
		nrm[n]=normalize(psi);
		for(int m=1;m<=expectSnaps;++m) //evolve in steps
		{
			if(timemethod==0 or timemethod==2) //choose from one of the evolution methods
				{
	        		psi=exactApplyMPO(expH,psi,args);
				}
			else
			if(timemethod==1 or timemethod==3)
				{
	        		fitApplyMPO(psi,expH,psi,args);
				}
			if(timemethod==2) psi=exactApplyMPO(expHb,psi,args);
			if(timemethod==3) fitApplyMPO(psi,expHb,psi,args);

		}

		for(int j = 1; j <= Ntotal; ++j)
		{
			psi.position(j);
			//densR[j-1][n] = (dag(prime(psi.A(j),Site))*sites.op(Nop,j)*psi.A(j)).real();
			densR[j-1][n] = getlocalExpectation(j, Nop,sites, psi);
			if(mpstype==3 or mpstype==6)
				{
				//densdnR[j-1][n] = (dag(prime(psi.A(j),Site))*sites.op(Nopdn,j)*psi.A(j)).real();
				//densdensR[j-1][n] = (dag(prime(psi.A(j),Site))*sites.op(NNop,j)*psi.A(j)).real();
				densdnR[j-1][n] = getlocalExpectation(j, Nopdn,sites, psi);

				}
			if(maxpsite>1 or (mpstype==3 or mpstype==6)){
				densdensR[j-1][n] = getlocalExpectation(j, NNop,sites, psi);
			}

		}

		//Calculate currents for time evolution
		current[0][n]=-2*(-currentcoupling*phase1*getcorrelatorOverlap(poscurrenta,poscurrentb,createOp,destroyOp,sites,psi,psi)).imag(); //Calculates Adag(a)*a(b),
		current[1][n]=-2*(-currentcoupling2*phase2*getcorrelatorOverlap(poscurrent2a,poscurrent2b,createOp,destroyOp,sites,psi,psi)).imag(); 
		current[2][n]=-2*(-currentcoupling3*phase3*getcorrelatorOverlap(poscurrent3a,poscurrent3b,createOp,destroyOp,sites,psi,psi)).imag(); 
			//println("Current t=0 :",current[j-1][0]);

		overlapinit[n]=pow(abs(overlapC(initpsi,psi)),2);

		bond[n]=maxM(psi);
    		comptime[n] = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
		placeholder[n]=0;

		
	//Save to file all snapshots taken so far
	
	//Save inside file as follows: 
	//time where snapshots where taken. In each line of length numberSnaps
	// current bond dimension
	// computational time
	//current 1 measured
	//current 2 measured
	//current 3 measured
	//Overlap with initial state
	//Not used
	//density for each site
	
	  if(stepsave==1){
		  datasetSave=dataset+"A"+ftos(timearray[currentnumberSnaps],dez);
		  datasetSave=ReplaceAll(datasetSave,".","_");
          }
	  std::ofstream myfile ("Density"+datasetSave+prefix+".dat");
	  if (myfile.is_open())
	  {
	    for(int count = 0; count <= currentnumberSnaps; count ++){
		myfile << timearray[count];
		if(count<currentnumberSnaps)
		{
			myfile <<",";
		}
	    }
	    myfile << "\n";

	    for(int count = 0; count <= currentnumberSnaps; count ++){
		myfile << bond[count];
		if(count<currentnumberSnaps)
		{
			myfile <<",";
		}
	    }
	    myfile << "\n";

	    for(int count = 0; count <= currentnumberSnaps; count ++){
		myfile << comptime[count];
		if(count<currentnumberSnaps)
		{
			myfile <<",";
		}
	    }
	    myfile << "\n";

	    for(int count = 0; count <= currentnumberSnaps; count ++){
		myfile << current[0][count];
		if(count<currentnumberSnaps)
		{
			myfile <<",";
		}
	    }
	    myfile << "\n";

	    for(int count = 0; count <= currentnumberSnaps; count ++){
		myfile << current[1][count];
		if(count<currentnumberSnaps)
		{
			myfile <<",";
		}
	    }
	    myfile << "\n";

	    for(int count = 0; count <= currentnumberSnaps; count ++){
		myfile << current[2][count];
		if(count<currentnumberSnaps)
		{
			myfile <<",";
		}
	    }
	    myfile << "\n";

	    for(int count = 0; count <= currentnumberSnaps; count ++){
		myfile << overlapinit[count];
		if(count<currentnumberSnaps)
		{
			myfile <<",";
		}
	    }
	    myfile << "\n";

	    for(int count = 0; count <= currentnumberSnaps; count ++){
		myfile << placeholder[count];
		if(count<currentnumberSnaps)
		{
			myfile <<",";
		}
	    }
	    myfile << "\n";
	    /*
	    for(int count = 0; count <= numberSnaps; count ++){
		myfile << entropy[count];
		if(count<numberSnaps)
		{
			myfile <<",";
		}
	    }
	    myfile << "\n";
	    */

	    for(int i = 0; i < Ntotal; i ++)
		{	
		    for(int count = 0; count <= currentnumberSnaps; count ++){
			myfile << densR[i][count];
			if(count<currentnumberSnaps)
			{
				myfile <<",";
			}
		    }
			myfile << "\n";
		}
	  if(mpstype==3 or mpstype==6) //For two species fermion down spin density
		{
		    for(int i = 0; i < Ntotal; i ++)
			{	
			    for(int count = 0; count <= currentnumberSnaps; count ++){
				myfile << densdnR[i][count];
				if(count<currentnumberSnaps)
				{
					myfile <<",";
				}
			    }
				myfile << "\n";
			}
		}
		//Density density correlations of spin up and down
	    for(int i = 0; i < Ntotal; i ++)
		{	
		    for(int count = 0; count <= currentnumberSnaps; count ++){
			myfile << densdensR[i][count];
			if(count<currentnumberSnaps)
			{
				myfile <<",";
			}
		    }
			myfile << "\n";
		}



	    myfile.close();
	    if(stepsave==1){
		    writeToFile("psi"+datasetSave,psi);
		    writeToFile("sites"+datasetSave,sites); //Write site object
		    if(currentnumberSnaps>1) cleanupoutput(datasetdelete,prefix);
		    datasetdelete=datasetSave;
	    }

	  }
	  else cout << "Unable to open file";

		if(mpstype==3 or mpstype==6)
			{
			printfln("Timestep= %.4f, Rtime=%.4f,Dsource=%.4f, Dsourcedn=%.4f, Ddrain=%.4f, Ddraindn=%.4f, overlap=%.4f,norm=%.4f, bond=%d",timearray[n],comptime[n],densR[0][n],densdnR[0][n],densR[Ntotal-1][n],densdnR[Ntotal-1][n],overlapinit[n],nrm[n],bond[n]);
			}
		else 
			{
			printfln("Timestep= %.4f, Rtime=%.4f, Dsource= %.4f, Ddrain=%.4f, overlap=%.4f, norm=%.4f, bond=%d",timearray[n],comptime[n],densR[0][n],densR[Ntotal-1][n],overlapinit[n],nrm[n],bond[n]);

			}

	}

    	
	  delete [] densR;
    }




    if(stepsave==0){
	    writeToFile("psi"+dataset,psi);
	    writeToFile("sites"+dataset,sites); //Write site object
    }
    duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;

    std::cout<<"Duration: "<< duration <<'\n';

    println("Program finished at time: ",finishtime);

    return 0;
    }
