///////////////////////////////////////////////////////////////////////////////
//
// THIS SOFTWARE IS PROVIDED "AS-IS". THERE IS NO WARRANTY OF ANY KIND.
// NEITHER THE AUTHORS NOR THE OHIO STATE UNIVERSITY WILL BE LIABLE
// FOR ANY DAMAGES OF ANY KIND, EVEN IF ADVISED OF SUCH POSSIBILITY.
//
// Copyright (c) 2010 Jyamiti Research Group.
// CS&E Department of the Ohio State University, Columbus, OH.
// All rights reserved.
//
// Author: Sayan Mandal
//
///////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <ctime>
#include <map>
#include <list>
#include <vector>
#include <cmath>
#include <sstream>
#include <unordered_set>
#include <string>
#include <sys/stat.h>

#ifndef INFINITY
#define INFINITY std::numeric_limits< double >::infinity()
#endif // INFINITY

#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
#include <boost/timer.hpp>
#include <boost/progress.hpp>

#include <CGAL/Cartesian.h>

#include<ParseCommand.h>
#include <Legal.h>
#include "SimplicialComplex.h"
#include "Complex.h"
#include <OFF_input_file.h>
#include <OFF_output_file.h>

// #ifndef COMPLEX_H
// #define COMPLEX_H

// #include "Complex.h"// header body...

// #endif

#include <Exception.h>


using namespace std;
using namespace trackLoop;
using namespace CGAL;
using namespace boost;
using namespace filesystem;
using namespace program_options;
using namespace Headers;

extern std::vector<int> complexSizes;
extern std::vector<int> accumulativeSizes;
extern int accumulativeSize;
extern float fThreshold;
extern vector<float> vecFiltrationScale;
extern SimplicialTree<bool> domain_complex;
//extern SimplicialTree<bool> range_complex;
//
//timer
extern std::clock_t start, timer1, timer2;
extern double dFuncTimeSum;
extern double dInsertTime;
extern double dCollapseTime;
extern int max_dimension;
extern int collapseCount;
extern int smallCount;

// this is the distance matrix used for matrix input
vector< vector<double> > g_distance_matrix;

// this is the adjacency list used for graph input
vector< vector< pair<int, double> > > g_adjacency_list;

// this is set to 0 for OFF input and 1 for MATRIX input.
int input_type_flag;

// template< typename Kernel_ > Complex< Kernel_ > complex, 

// outer vector: loop, inner vector: edges forming the loop
typedef std::vector<std::vector<int>> higherOrder;
typedef std::vector<std::vector<std::vector<float>>> higherPoint;

// void createTreeLoopTracker_ptr(std::vector<std::vector<std::vector<int>>> simplex_vertices);
bool CheckBoundary (std::vector<std::vector<int>> simplex_vertices);
void CheckBoundaryBirthOfLoops(std::map<int, higherOrder> birthOfLoops);
// ListNodeptr AddBoundary(std::vector<std::vector<int>> simplex_vertices);
// int independantCycleCalculate2(std::map<float,int> weights, std::vector<std::vector<std::vector<int>>> higher_simplex);
int independantCycleCalculate(std::multimap<float, higherOrder> &simp_weight);


int lowfromHigherOrder(higherOrder &simplex_vertices);

void ComputingPersistenceForSimplicialMapElementary(vector<vector<int> > simplex, const char* file_name_of_domain_complex,
	bool is_domain_complex_with_annotation,
	vector<string>& vecElemOpers,
	bool is_save_range_complex_with_annotation = false,
	const char* new_range_complex_file_name = NULL);

bool bornTracker(higherOrder simp2, std::map<int, higherOrder> birthOfLoops);


bool barcodeCompare(const pair<int, int>& a, const pair<int, int>& b)
{
	float a1, a2, b1, b2;
	a1 = vecFiltrationScale[a.first];
	a2 = a.second == -1 ? vecFiltrationScale.back() : vecFiltrationScale[a.second];
	b1 = vecFiltrationScale[b.first];
	b2 = b.second == -1 ? vecFiltrationScale.back() : vecFiltrationScale[b.second];
	return a2 - a1 < b2 - b1;

}

bool barcodeCompareWithTS(const std::pair<int, pair<int, int>>& a, const std::pair<int, pair<int, int>>& b)
{
	return barcodeCompare(a.second, b.second);
}


void printHigherOrder(higherOrder ho){
	// typedef std::vector<std::vector<int>> higherOrder;
	for(int i=0;i<ho.size();i++){
		for(int j=0;j<ho[i].size();j++)
			cout<<ho[i][j]<<" ";
		cout<<"<-->";
	}
	cout<<"\n";
}


void printMultimap(std::multimap<float, higherOrder> simp_weight, int op){
	int countitsw = 0;
	for (std::multimap<float, higherOrder>::iterator itsw = simp_weight.begin(); itsw != simp_weight.end(); itsw++){
		if (countitsw == op)
		{
			printHigherOrder(itsw->second);
			return;
		}
		countitsw++;
	}

}

bool checkPoints(std::vector<float> a, std::vector<float> b){
	for(int i=0;i<a.size();i++)
		if(a[i]!=b[i])
			return false;
	return true;
}

higherPoint modifylastloopPoint(higherPoint lastLoop){
	higherPoint temp;
	temp.push_back(lastLoop[0]);
	lastLoop.erase(lastLoop.begin());

	int count = 1;
	int fullsize = lastLoop.size();
	while(count<=fullsize){
		std::vector<float> a;
		a.push_back(temp[temp.size()-1][1][0]);
		a.push_back(temp[temp.size()-1][1][1]);
		a.push_back(temp[temp.size()-1][1][2]);
		// Point a = temp[temp.size()-1][1];	//last element in last edge of vector, last vertex;
		int k;
		for(k=0;k<lastLoop.size();k++)	{
			if(checkPoints(lastLoop[k][1],a)==true ||checkPoints(lastLoop[k][0],a)==true)
				break;
		}
		if(checkPoints(lastLoop[k][1],a))		//need to reverse the edge
			{	std::vector<std::vector<float>> tttt;
				tttt.push_back(lastLoop[k][1]);
				tttt.push_back(lastLoop[k][0]);
				temp.push_back(tttt);
			}
		else
			temp.push_back(lastLoop[k]);
		lastLoop.erase(lastLoop.begin() + k);
		// cout<<"2";
		count++;
	// // getchar();

	}
	return temp;
	
}

higherOrder modifylastloop(higherOrder lastLoop){
	higherOrder temp;
	for(int i=0;i<lastLoop.size();i++){
		// for(int j=0;j<lastLoop[i].size();j++){
			if(lastLoop[i].size()!=2)
				{cout<<"error"; exit(0);}
			cout<<lastLoop[i][0]<<" "<<lastLoop[i][1];
			cout<<"-+-";
		}

	temp.push_back(lastLoop[0]);
	lastLoop.erase(lastLoop.begin());
	int count = 1;
	int fullsize = lastLoop.size();
	while(count<=fullsize){
		int a = temp[temp.size()-1][1];	//last element in last edge of vector, last vertex;
		int k;
		for(k=0;k<lastLoop.size();k++)
			if(lastLoop[k][1]==a ||lastLoop[k][0]==a)
				break;
			else if(k==lastLoop.size()-1)
			{
				cout<<"Loop mismatch: "<<a<<" ";
				for (int it = 0; it < temp.size(); ++it)
					cout<<temp[it][0]<<" "<<temp[it][1]<<"|";
				cout<<"\n";
				for (int it = 0; it < lastLoop.size(); ++it)
					cout<<lastLoop[it][0]<<" "<<lastLoop[it][1]<<"|";
				exit(0);
			}

		if(lastLoop[k][1]==a)		//need to reverse the edge
			{	
				std::vector<int> tttt;
				tttt.push_back(lastLoop[k][1]);
				tttt.push_back(lastLoop[k][0]);
				temp.push_back(tttt);
				
			}
		else if(lastLoop[k][0]==a){
			cout<<"higherOrder size 5:"<<lastLoop.size();
			// getchar(); 
			temp.push_back(lastLoop[k]);
		}
		else
		{
			cout<<"Loop mismatch"; exit(0);
		}

			cout<<"higherOrder size 3:"<<lastLoop.size();
	// getchar();

		lastLoop.erase(lastLoop.begin() + k);
		// cout<<"2";
		count++;
	// getchar();

	}
	return temp;
	
}

void loopPrinting( std::map<int,higherPoint> vloop_all, std::string loops_folder){

//higherpoint: vector, vector , vector float
	boost::filesystem::path dir(loops_folder.c_str());
	boost::filesystem::create_directory(dir);
	// std::vector<higherOrder> vInd;
	// const int dir_err = mkdir(loops_folder.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
	// if (-1 == dir_err)
	// {
	//     printf("Error creating directory!n");
	//     exit(1);
	// }
	// chmod("./myfile", S_IRWXU);
	higherOrder vInd;
	
	// for(int l1=0;l1<vloop.size();l1++)
	for(std::map<int,higherPoint>::iterator iter = vloop_all.begin(); iter != vloop_all.end(); ++iter)
	{
		std::vector<int> edge;
		int k =  iter->first;
		higherPoint vloop = iter->second;
	// {
		std::string file = loops_folder+std::to_string(k)+".off";
		cout<<"OFF folder:"<<file<<"\n";
		ofstream ofloop(file.c_str());
		ofloop<<"OFF"<<std::endl<<vloop.size()*4<<" "<<vloop.size()<<" 0\n";
		int count = 0;
		for(int l2=0;l2<vloop.size();l2++){
			// if(vloop[l1][l2][0].size()!=3 || vloop[l1][l2][1].size()!=3){
			// cout<<"Cannot print loop, not 3D";
			// exit(0);
			// }
		ofloop<<vloop[l2][0][0]<<" "<<vloop[l2][0][1]<<" "<<vloop[l2][0][2]<<"\n";
		ofloop<<vloop[l2][1][0]<<" "<<vloop[l2][1][1]<<" "<<vloop[l2][1][2]<<"\n";
		ofloop<<vloop[l2][0][0]<<" "<<vloop[l2][0][1]<<" "<<vloop[l2][0][2]<<"\n";
		ofloop<<vloop[l2][1][0]<<" "<<vloop[l2][1][1]<<" "<<vloop[l2][1][2]<<"\n";

		edge.push_back(count++);
		edge.push_back(count++);
		edge.push_back(count++);
		edge.push_back(count++);
		vInd.push_back(edge);
		edge.clear();
		// ofloop<<"#######################################################\n";

		}

		for(int ed=0;ed<vInd.size();ed++)
		{
			ofloop<<"4 "<<vInd[ed][0]<<" "<<vInd[ed][1]<<" "<<vInd[ed][2]<<" "<<vInd[ed][3]<<" 1.0 1.0 0.0"<<std::endl;
		}

		ofloop.close();
		vInd.clear();

	}
	// count = 0;

	
}


int simpersPart(std::vector<int>  &born,std::vector<int>  &dead, std::string simpers_file){

	ifstream ff(simpers_file.c_str());
    // cout<<"simpers part";
    if( ff.good()==false)
    {
        cout<<"simpers file "<<simpers_file<<" does not exist.";
        exit(0);
    }

    int dim, fborn;
    std::string fdead;
    while(!ff.eof()){

    	char sLine[256]="";
    	ff.getline(sLine, 256);
    	if(sLine==""||strlen(sLine)==0)
    		return 0;
    	// cout<<sLine<<",,,, ";
    	stringstream ss;
		ss.str(sLine);
			
    	ss >> dim;	ss>>fborn; ss>> fdead;
    	
    	if(dim==0)
    		continue;
    	if(dim==2)
    		return 0;
    	born.push_back(fborn);
    	// cout<<fborn<<"fb ";
    	if(fdead == "inf")
    		dead.push_back(-1);
    	else
    			dead.push_back(stoi(fdead));

    	// cout<<fdead<<"fd"<<"\n";

    }
    return 0;
}


int loopExistenceChecker(higherOrder simp2, std::map<int, higherOrder> birthOfLoops){
	//true: it already exists. danger
	// false: new one does not exist.
	// int iterator=0;
	for(std::map<int,higherOrder>::iterator iter = birthOfLoops.begin(); iter != birthOfLoops.end(); ++iter)
	{
		bool flag = false;
		higherOrder buff = iter->second ;
		if(buff.size()!=simp2.size())
			continue;
		for(int ord=0;ord<buff.size();ord++){
			if(buff[ord][0]!=simp2[ord][0] || buff[ord][0] != simp2[ord][1] )
				break;
			else if(ord==buff.size()-1)
				return iter->first;
		}
		// iterator++;

	}
	return -1;
}// end of loopExistenceChecker

int main( int argc, char *argv[] )
{
	
	std::string input_pointcloud_file;
	// std::string output_persistence_file_name;
	std::string loops_folder;
	std::string filtration_file;
	// std::string output_file;
	std::string point_file;
	std::string simpers_file;
	string input_file_name;
	double alpha=INFINITY;
	bool bDeathTimeOrder = true;
	bool bTimeStamp = false;
	float fMaxScale = 0.0;
	double sampling_coefficient = 1;
	// float born, dead;
	std::vector<int> vborn;
	std::vector<int> vdead;
	std::map<int, higherPoint> vloop;	//Point, two points: one edge, set of edge: one cycle, set of cycle: against birth-time

	int dimensions, noPoints; 
	float scalecount = 0;
	std::vector<Point> allPts;
	higherOrder storeComplexforShortLoop; // stores entire complex to be pushed into shortloop
	std::map<int, higherOrder> birthOfLoops;	//int: birth time, higherOrder: edges in the loop
	std::map<int, int> nedges; //number of edges


	ParseCommand(argc, argv, input_pointcloud_file, 
		filtration_file, sampling_coefficient);
	// born = 100;
	// dead = 100;

    // ifstream f(input_pointcloud_file.c_str());
    simpers_file = input_pointcloud_file.substr(0,input_pointcloud_file.size()-4)+"pers.txt"; //Input to simpers
    loops_folder = input_pointcloud_file.substr(0,input_pointcloud_file.size()-4)+"loops/"; //Output of loops
    ifstream pf(input_pointcloud_file.c_str());

    
    
    if( pf.good()==false)
    {
        cout<<"Point file does not exist.";
        return 0;
    }

	ifstream ff(filtration_file.c_str());
	// string barnumber = filtration_file+"_pers";
	// ofstream ff2((filtration_file+"_pers").c_str());
	if(ff.good()==false)
    {
        cout<<"Filtration file does not exist.";
        exit(0);
    }

    typedef Cartesian< double > Kernel;
    

    simpersPart(vborn, vdead, simpers_file);
    // vloop.reserve(vborn.size());
	//OFF_input_file< Kernel > input( input_pointcloud_file.c_str() );    
	// vborn, vdead, 
	

	pf >> dimensions ;
	pf >> noPoints;
	cout<<"Dim: "<<dimensions<<" #Pt: "<<noPoints<<" \n";
	bool verbose = true;
	// Complex< Kernel > complex( dimensions, verbose );		
	const char* buff = '\0';

	vector<string> vecElemOpers;
	vector<int> interm;
	
	
	// Create vertices SHORTLOOP

	domain_complex.bGenerator = false;
	for ( int itp=0; itp < noPoints; itp++ )
	{			
		Point p(dimensions);
		for(int k=0;k<dimensions;k++){
			float coord;
			pf >> coord;
			p.set_coord(k,coord);
		}
		
		// Vertex< Kernel > *p_vertex( new Vertex< Kernel >(p));
		// complex.insert_vertex( p_vertex );
		scalecount+=1;
		filtration_step += 1;
		timer2 = std::clock();
		allPts.push_back(p);
		vecFiltrationScale.push_back(scalecount);
		interm.push_back(itp);

		domain_complex.ElementaryInsersion(interm);
		interm.clear();
		dFuncTimeSum += (std::clock() - timer2);

		complexSizes.push_back(domain_complex.ComplexSize());
		accumulativeSizes.push_back(domain_complex.accumulativeSimplexSize);

	}

	// Add edges and triangles; edges are created if missing
	
	bool bornflag = false;
	bool deadflag = false;
	int currentv1=-1,currentv2=-1;
    while (!ff.eof())	
    {

		// if(deadflag==true){
			// cout<<"Got here dead";
			// getchar();
		// }

    	char sLine[256]="";

    	ff.getline(sLine, 256);
    	// cout<<"sLine:"<<sLine<<"\n";
    	if(sLine[0]=='c'||strlen(sLine)==0)
    		continue;


    	stringstream ss;
		ss.str(sLine);
		char ic;
		int index;
		float indf;
		std::vector<int> simplex1;

		ss >> ic;
		
		if(ic=='#'){
			ss >> indf;

			auto itf = std::find(vdead.begin(), vdead.end(), indf);
			// ******************** DEAD PART *********************
			if(itf != vdead.end() ){	
				cout<<"Short Loop Dead:sL:"<<sLine<<"\n"; 
				// getchar();
				// CheckBoundaryBirthOfLoops(birthOfLoops);
				int index = std::distance(vdead.begin(), itf);
				int lid = vborn[index]; //loop_index_which_died
				higherOrder lwd = birthOfLoops[lid]; // actual loop which died
				// cout<<"dead: "<<lid;
				// getchar();
				if(nedges[lid]!=lwd.size()){
					cout<<"Number mismatch of loops: "<<nedges[lid]<<" "<<lwd.size();
					exit(0);
				}
				lwd = modifylastloop(lwd);
				cout<<"Loop born at: "<<vborn[index]<<", died at: "<<vdead[index]<<"\n";
				// cout<<"dead 2";
				// getchar();
				printHigherOrder(lwd);
				// getchar();
				// cout<<"erase";
				birthOfLoops.erase(lid);
				nedges.erase(lid);
				deadflag = true;
				cout<<"short loop second";
				CheckBoundaryBirthOfLoops(birthOfLoops);
				// getchar();
				continue;
			}
			// ******************* BORN PART ************************
			else if (std::find(vborn.begin(), vborn.end(), indf) != vborn.end()){
				// {cout<<"Neither constructor or destructor."<<std::endl;  continue;}	// else{
				
				
				cout<<"Short Loop Born: "<<indf<<"|simplex: ";
				CheckBoundaryBirthOfLoops(birthOfLoops);
				bool loopadded = false;
				Complex< Kernel > complex( dimensions, verbose );
    	
		    	for ( int itp=0; itp < noPoints; itp++ )
				{			
					Point p(dimensions);
					p = allPts[itp];
					Vertex< Kernel > *p_vertex( new Vertex< Kernel >(p));
					complex.insert_vertex( p_vertex );
					
				}
				for( int itp=0;itp<storeComplexforShortLoop.size();itp++){
					if(storeComplexforShortLoop[itp].size()==2){
						Vertex< Kernel > &a( complex.vertex_at( storeComplexforShortLoop[itp][0] ) );
						Vertex< Kernel > &b( complex.vertex_at( storeComplexforShortLoop[itp][1] ) );
						// cout<<a.index()<<"... "<<b.index()<<"|";
						complex.create_edge( a, b, sqrt(a.location().get_squared_distance_to(b.location())) );
					}
					else if(storeComplexforShortLoop[itp].size()==3){
						Vertex< Kernel > &a( complex.vertex_at( storeComplexforShortLoop[itp][0] ) );
						Vertex< Kernel > &b( complex.vertex_at( storeComplexforShortLoop[itp][1] ) );
						Vertex< Kernel > &c( complex.vertex_at( storeComplexforShortLoop[itp][2] ) );
						// cout<<a.index()<<"..."<<b.index()<<"..."<<c.index()<<"|";
						complex.create_triangle( a, b, c );
					} 	

				}
				cout<<"\n";

				cout << complex.number_of_vertices() << " vertices" << endl;
				cout << complex.number_of_edges() << " edges" << endl;
				cout << complex.number_of_triangles() << " triangles" << endl;
				cout << complex.number_of_vertices() + complex.number_of_edges()
				+ complex.number_of_triangles() << " simplices total" << endl;
				complex.contract();	// builds the tree
				complex.sample( sampling_coefficient );		// Gets a random sample from the complex.  All points are used if sampling_coefficient=1.
				
				complex.compute_basis();
				
				cout << complex.basis_rank() << " loops\n";
				//Find which loop is born here
				higherOrder simp2;
				higherPoint vP;
				cout<<"currents: "<<currentv1<<" "<<currentv2<<std::endl;
				for ( unsigned i( 0 ); i != complex.basis_rank(); ++i )
				{

					// if(loopadded==true && indf!=564)
						// break;
					Basis_loop< Kernel > &basis_loop = complex.basis_loop_at( i ) ;
					cout << "Loop " << i << " (" << basis_loop.size();
					cout << " edges, length=" << basis_loop.norm() << "):";
					Basis_loop< Kernel >::Iterator it_edge( basis_loop.begin() );
					std::vector<int> interm;
					std::vector<std::vector<float>> vP2;
					// getchar();
					bool containsthisedge = false;
					for ( ; it_edge != basis_loop.end(); ++it_edge )
					{
						Edge< Kernel > &edge( **it_edge );
						cout<<edge.a().index()<<" "<<edge.b().index()<<"...";
						if((currentv1==edge.a().index() && currentv2==edge.b().index())||(currentv2==edge.a().index() && currentv1==edge.b().index()))
							containsthisedge = true;
						// else containsthisedge = false;
						interm.push_back(edge.a().index());
						interm.push_back(edge.b().index());
						simp2.push_back(interm);
						interm.clear();

						int dimhere = edge.a().m_location.get_dim();

						std::vector<float> coord;
						for(int idim=0;idim<dimhere;idim++)
							coord.push_back(edge.a().m_location.get_coord(idim));
						
						vP2.push_back(coord);
						
						coord.clear();
						for(int idim=0;idim<dimhere;idim++)
							coord.push_back(edge.b().m_location.get_coord(idim));

						vP2.push_back(coord);
						
						vP.push_back(vP2);
						vP2.clear();
					}
					// getchar();
					cout<<"\n";

					if(simp2.size()!=basis_loop.size()){
						cout<<"Not all edges are inserted: "<<simp2.size()<<" "<<basis_loop.size()<<"\n";
						printHigherOrder(simp2);
						exit(0);
					}

					// Takes in the current loop and set containing all loops and 
					// sees if this current one is independant, if so then this was born
					if( containsthisedge==true && bornTracker(simp2,birthOfLoops)==true){
						// cout<<"no it aint";
						// getchar();
						loopadded = true;
						// if(indf==564)
							// {cout<<"chosen: ";printHigherOrder(simp2); }
						// getchar();}
						int indy = loopExistenceChecker(simp2, birthOfLoops);
						if(indy!=-1)
						{

							cout<<"Wrong born tracker: "<<std::endl;
							printHigherOrder(simp2);
							cout<<"**********************"<<"\n";
							printHigherOrder(birthOfLoops[indy]);
							exit(0);
						}
						
						birthOfLoops.insert(std::pair<int, higherOrder>(indf, simp2));
						nedges.insert(std::pair<int, int>(indf, basis_loop.size()));
						if(simp2.size()!=basis_loop.size()){
							cout<<"Number mismatch in loops even before insert";
							exit(0);
						}

						vloop.insert(std::pair<int, higherPoint>(indf, vP));
					}
						simp2.clear();
						vP.clear();
		
				}//"Loop basis rank"
				// cout<<"out of basis part";
				if(loopadded == false)
				{
					cout<<"No loop has been added. This is an error";
					exit(0);
				}
				
				cout<<"born end";
				CheckBoundaryBirthOfLoops(birthOfLoops);
				continue;
			}// Born Part

			else continue;// Neither constructor nor destructor
		}// if ic = # part	
		// if(deadflag==true)
		// {
		// 	cout<<"got here";
		// 	getchar();
		// }

		while (ss >> index)
		{
		simplex1.push_back(index);
		}

		if(simplex1.size() <= 1)
		{
		cout<<" vertices should already be inside. Error\n";
		exit(0);
		}
		
		if(simplex1.size()==2)
		{
			currentv1 = simplex1[0];
			currentv2 = simplex1[1];
		}

		storeComplexforShortLoop.push_back(simplex1);
		scalecount+=1;
		filtration_step += 1;
		timer2 = std::clock();

		vecFiltrationScale.push_back(scalecount);
		domain_complex.ElementaryInsersion(simplex1);
		simplex1.clear();
			
		complexSizes.push_back(domain_complex.ComplexSize());
		accumulativeSizes.push_back(domain_complex.accumulativeSimplexSize);
		dFuncTimeSum += (std::clock() - timer2);
		
    }//end of while going through each line of the simpers code

// cout<<"this part";
// getchar();

loopPrinting(vloop, loops_folder);
return 1;

}// end of main