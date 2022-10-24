#include "itensor/all.h"
#include "itensor/util/print_macro.h"
#include <fstream>
#include <iostream>
using namespace itensor;
using std::vector;

/*
void collectdata(float t, float sy, float sx, float overlp)
{
    std::fstream outfile("out_L=81_D=0_theta=30.txt", std::ios::app | std::ios::in);
    //	std::fstream outfile("try.txt",std::ios::app | std::ios::in);
    outfile << t << " "<< sy <<" "<< sx << " " << overlp << "\n";
    outfile.close();
}
*/
int main(int argc, char *argv[])
{
    int N = atoi(argv[1]);
    double D = atof(argv[2]);
    double theta = atof(argv[3]); //in degree
    double pi = 2 * acos(0.0);
    theta = theta*pi/180;

    printfln("N= ", N);

    double Jf = 1.0;
    //double D = 0.35;

    printfln("Jf= ",Jf);
    printfln("D= ",D);
    printfln("theta (in rad) = ",theta);
   
    auto sites = SpinOne(N, {"ConserveQNs=", false});

    //rotating the state with acting field at theta angle
    auto ampo = AutoMPO(sites);
    for (int j = 1; j <= N; j++)
    {
        ampo += -cos(theta), "Sz", j;
        ampo += -sin(theta), "Sy", j;
    }
    auto H = toMPO(ampo);

    //set initial wavefunction matrix product state
    auto state = InitState(sites);

    for (int i = 1; i <= N; i++)
    {
        //if (i<=N/2){state.set(i, "Up");}
	//else {state.set(i, "Dn");}
	state.set(i, "Up");
    } //initial state
    auto phi = MPS(state);

    //auto phi = randomMPS(sites);
    
    // Sweeping of DMRG
    auto sweeps = Sweeps(100);
    sweeps.maxm() = 5, 10,20, 50, 100, 100, 200, 200, 200, 500, 500, 1000;
    //sweeps.maxm() = 100,200,200,500,500,500,1000,1000,2000,2000,2000,3000,3000,3000;
    sweeps.cutoff() = 1E-12;
    sweeps.niter() = 2;
    sweeps.noise() = 1E-6, 1E-6, 1E-7, 1E-8, 0, 0, 0, 0;
    //println(sweeps);

    //Begin DMRG calculation
    auto [energy, psi] = dmrg(H, phi, sweeps, "Quiet");
  
    printfln("norm= ", norm(psi));
    printfln("initial FM abs(overlapC(psi,phi)))= ",abs(overlapC(psi,phi)));
    printfln("Lowest energy = ", energy);

    Real tstep = 0.02; //time step (smaller is generally more accurate)
    //Real ttotal = 1.0; //total time to evolve
    Real cutoff = 1E-8; //truncation error cutoff when restoring MPS form

    //Create a std::vector (dynamically sizeable array)
    //to hold the Trotter gates
    auto gates = vector<BondGate>();

/*
#ifdef EVEN_ODD
    for (int b = 1 ; b <= N-1 ; b += 2)
    {
        auto hterm = -Jf * op(sites, "Sz", b) * op(sites, "Sz", b + 1);
        hterm += -0.5 * Jf * op(sites, "S+", b) * op(sites, "S-", b + 1);
        hterm += -0.5 * Jf * op(sites, "S-", b) * op(sites, "S+", b + 1);

        auto g = BondGate(sites, b, b + 1, BondGate::tReal, tstep / 2., hterm);
        gates.push_back(g);

    }
#endif
*/
    //Create the gates exp(-i*tstep/2*hterm)
    //and add them to gates
    for (int b = 1; b <= N - 1; ++b)
    {
        auto hterm = -Jf * op(sites, "Sz", b) * op(sites, "Sz", b + 1);
        hterm += -0.5 * Jf * op(sites, "S+", b) * op(sites, "S-", b + 1);
        hterm += -0.5 * Jf * op(sites, "S-", b) * op(sites, "S+", b + 1);

        auto g = BondGate(sites, b, b + 1, BondGate::tReal, tstep / 2., hterm);
        gates.push_back(g);
    }

    //Create the gates exp(-i*tstep/2*hterm) in reverse
    //order (to get a second order Trotter breakup which
    //does a time step of "tstep") and add them to gates
    for (int b = N - 1; b >= 1; --b)
    {
        auto hterm = -Jf * op(sites, "Sz", b) * op(sites, "Sz", b + 1);
        hterm += -0.5 * Jf * op(sites, "S+", b) * op(sites, "S-", b + 1);
        hterm += -0.5 * Jf * op(sites, "S-", b) * op(sites, "S+", b + 1);

        auto g = BondGate(sites, b, b + 1, BondGate::tReal, tstep / 2., hterm);
        gates.push_back(g);
    }

    //for onsite terms use the identity operator to make it two site operator
    for (int b = 1; b <= N - 1; ++b)
    {
        auto hterm = -D * op(sites, "Sz2", b) * op(sites, "Id", b + 1);
        //hterm += -B*op(sites,"Sx",b)*op(sites,"Id",b+1);
        hterm += -D * op(sites, "Sz2", b + 1) * op(sites, "Id", b);
        //hterm += -B*op(sites,"Sx",b+1)*op(sites,"Id",b);

        if (b == 1) //for boundary,as for boundary, gates will be applied only once per site
        {
            hterm += -D * op(sites, "Id", b + 1) * op(sites, "Sz2", b);
            //hterm += -B*op(sites,"Id",b+1)*op(sites,"Sx",b);
        }
        if (b == N - 1)
        {
            hterm += -D * op(sites, "Id", b) * op(sites, "Sz2", b + 1);
            //hterm += -B*op(sites,"Id",b)*op(sites,"Sx",b+1);
        }
        auto g = BondGate(sites, b, b + 1, BondGate::tReal, tstep / 2., hterm);
        gates.push_back(g);
    }

    //Save initial state;
    auto psi0 = psi;
    auto psi0dag = dag(psi0);
    //auto sp_psi = op(sites,"S+",(N+1)/2)*psi((N+1)/2);
    //sp_psi.noPrime();
    //psi.set((N+1)/2,sp_psi);

    //Time evolve, overwriting psi when done
    printfln("time, <sx>, <sy>, <psi|psi(t)>");

    for (int t = 1; t <= 10000; ++t)
    {
        gateTEvol(gates, tstep, tstep, psi, {"Cutoff=", cutoff, "Verbose=", false});

        double sxt = 0;
	double syt = 0;

        for (int i = 1; i <= N; ++i)
        {
            psi.position(i);
            auto sx = op(sites, "Sx", i);
	    auto sy = op(sites, "Sy",i);

            auto Cx = psi.A(i);
	    auto Cy = psi.A(i);

            Cx *= sx * dag(prime(psi.A(i), "Site"));
	    Cy *= sy * dag(prime(psi.A(i), "Site"));

            sxt += eltC(Cx).real();
	    syt += eltC(Cy).real();
	    
            //printfln("",t, " ", eltC(C));
        }

        printfln("", t * tstep, " ", sxt, " ", syt," ", abs(innerC(psi0, psi)));
        //collectdata(t*tstep,sxt,syt,abs(innerC(psi0, psi)));
    }

    /*
if (i==1){C *= sm*prime(psi0dag(i),"Site");}
else {C *= psi0dag(1);}

for (int k=2;k<=N;++k)
{
if (k==i){C *= psi(k);
	C *= sm*prime(psi0dag(k),"Site");}
else {C *= psi(k)*psi0dag(k);}
}
}
//auto result = -Cplx_i*sqrt(2)*exp(Cplx_i*energy*t*tstep)*eltC(C);
//collectdata(t*tstep,i,result.real(),result.imag());

printfln("i= ",i," t= ",t*tstep," C= ", result);
//time.push_back(float(t*tstep));
//pos.push_back(float(i));
//fxt.push_back(float(result.imag()));
}
printfln("Maximum MPS bond dimension after time evolution is %d",maxLinkDim(psi));
}
printfln("Lowest energy = ",energy);
*/
    return 0;
}
