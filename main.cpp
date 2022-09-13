#include <vector>
#include <stdexcept>
#include <iostream>
#include <fstream>
#include "optimization.h"
#include "DeformationGraph.h"



// ongoing work for an obj file parser
std::vector<std::vector<std::vector<double> > > parseFile()
{
	std::vector<std::vector<double> > coord;
	std::vector<std::vector<double> > coord_constraint;
	std::vector<std::vector<std::vector<double> > > data;
	std::ifstream reader;
	//~ reader.open("cube.obj");
	reader.open("canonical-shape.obj");
	std::vector<double> x,y,z;
	if (reader.is_open()) 
	
	{
		std::string line;
		while (std::getline(reader, line))
		{
			if (line.substr(0,2) != "v ")
				continue;
			line = line.substr(2);
			std::string::size_type sz;     // alias of size_t
			x.push_back(std::stod(line, &sz));
			line = line.substr(sz);
			y.push_back(std::stod(line, &sz));
			line = line.substr(sz);
			z.push_back(std::stod(line, &sz));
		}
		//~ for (int i=0;i<x.size();i++){
			//~ std::cout << "x:" << x.at(i) << ", y:" << y.at(i) << ", z:" << z.at(i) << std::endl;
		//~ }
	}
	else {
		std::cout << "could not open input.obj" << "\n";
	}
	coord.push_back(x);
	coord.push_back(y);
	coord.push_back(z);
	reader.close(); 
	data.push_back(coord);
	
	//~ reader.open("cube_modif.obj");
	//~ reader.open("kinect-capture1.obj");
	reader.open("canonical-shape_test.obj");
	//~ reader.open("canonical-modif.obj");
	std::vector<double> x1,y1,z1;
	if (reader.is_open()) 
	
	{
		std::string line;
		while (std::getline(reader, line))
		{
			if (line.substr(0,2) != "v ")
				continue;
			line = line.substr(2);
			std::string::size_type sz;     // alias of size_t
			x1.push_back(std::stod(line, &sz));
			line = line.substr(sz);
			y1.push_back(std::stod(line, &sz));
			line = line.substr(sz);
			z1.push_back(std::stod(line, &sz));
		}
		//~ for (int i=0;i<x.size();i++){
			//~ std::cout << "x:" << x.at(i) << ", y:" << y.at(i) << ", z:" << z.at(i) << std::endl;
		//~ }
	}
		else {
			std::cout << "could not open input.obj" << "\n";
		}
		coord_constraint.push_back(x1);
		coord_constraint.push_back(y1);
		coord_constraint.push_back(z1);
		reader.close(); 
		data.push_back(coord_constraint);
		
		return data;
}


int main(int , char** )
{
	int nb_nodes=20;
	int k_nearest=4;

	//~ int nb_nodes=10;
	//~ int k_nearest=2;	
	
	std::vector<std::vector<std::vector<double> > > data;
	
	data=parseFile();
	std::vector<std::vector<double> > coord,coord_constraint;
	coord=data[0];
	coord_constraint=data[1];
	int nb_vertices=coord.at(0).size();	
	int nb_vertices_modif=coord_constraint.at(0).size();
	
	//~ std::cout<<nb_vertices<<"\n";
	//~ std::cout<<nb_vertices_modif<<"\n";
	//~ for (int i=0; i<nb_vertices;i++){
	//~ std::cout << "x:" << coord.at(0).at(i)<< ", y:" << coord.at(1).at(i) << ", z:" << coord.at(2).at(i) << std::endl;
	//~ }
	
    double* v[nb_vertices];
    
    for (int i = 0 ; i < nb_vertices ; i++)
    {
		
        v[i] = new double[3];
        for (int j = 0 ; j < 3 ; j ++) v[i][j] = coord.at(j).at(i);
    }

	double** human = v;

	std::vector<std::vector<double> >nodes;
	std::vector<double> x_nodes,y_nodes,z_nodes;
	int k=0;
	int step=nb_vertices/(nb_nodes);
	//~ std::cout<<"step= "<<step<<std::endl;
	
	for (int i = step-1 ; i < nb_vertices ; i+=step)
	{		     
        x_nodes.push_back(v[i][0]);
        y_nodes.push_back(v[i][1]);
        z_nodes.push_back(v[i][2]);
	}
	nodes.push_back(x_nodes);
	nodes.push_back(y_nodes);
	nodes.push_back(z_nodes);
        		
    
    int count=0;
    double eps=1e-3;
    std::vector<int> indice;
	for (int i = 0 ; i < nb_vertices_modif ; i++)
	{	
        if ((coord_constraint.at(0).at(i)-coord.at(0).at(i)>eps) ||(coord_constraint.at(0).at(i)-coord.at(0).at(i)<-eps))
        {
			std::cout<<"x_cont!=x0"<<std::endl;
			indice.push_back(i);
			count++;
		}
		else if ((coord_constraint.at(1).at(i)-coord.at(1).at(i)>eps) ||(coord_constraint.at(1).at(i)-coord.at(1).at(i)<-eps))
		{
			std::cout<<"y_cont!=y0"<<std::endl;
			indice.push_back(i);
			count++;
		}
		else if ((coord_constraint.at(2).at(i)-coord.at(2).at(i)>eps) ||(coord_constraint.at(2).at(i)-coord.at(2).at(i)<-eps))
		{
			std::cout<<"z_cont!=z0"<<std::endl;
			indice.push_back(i);
			count++;
		}
		else if(coord_constraint.at(2).at(i)>=0)
		{
			//~ std::cout<<"z_cont<0"<<std::endl;
		}
		else 
		{
			//~ std::cout<<i<<std::endl;
		}
			 
			 
			
	}
    
    
	//~ for (int i = step-1 ; i < nb_vertices ; i+=step)
	//~ {	
		//~ std::cout<<"indice i : "<<i<<" ----- "; 	
        //~ for (int j = 0 ; j < 3 ; j ++) std::cout<<v[i][j]<<" ";
		//~ std::cout<<"\n";        
    //~ }
   	std::cout<<"count= "<<count<<"\n";
	//~ for (int i=0; i<nb_vertices_modif;i++){
	//~ std::cout << "x:" << coord_constraint.at(0).at(i)<< ", y:" << coord_constraint.at(1).at(i) << ", z:" << coord_constraint.at(2).at(i) << std::endl;
	//~ }
	
	//~ for (int i=0; i<count;i++){
	//~ std::cout << "indice:" << indice.at(i)<< std::endl;
	//~ }
	int fixed=10;
    double* q[count+fixed];
    
    for (int i = 0 ; i < count ; i++) 
    {
		double* q0 = new double[3];
		
		q0[0]=coord_constraint.at(0).at(indice.at(i));
		q0[1]=coord_constraint.at(1).at(indice.at(i));
		q0[2]=coord_constraint.at(2).at(indice.at(i));

		q[i]=q0;

	}
	for (int i = 0 ; i < 5 ; i++) 
    {
		double* qf = new double[3];
		qf[0]=coord_constraint.at(0).at(4032+i);
		qf[1]=coord_constraint.at(1).at(4032+i);
		qf[2]=coord_constraint.at(2).at(4032+i);
		q[count+i-1]=qf;
	}
		for (int i = 0 ; i < 5 ; i++) 
    {
		double* qf = new double[3];
		qf[0]=coord_constraint.at(0).at(1568+i);
		qf[1]=coord_constraint.at(1).at(1568+i);
		qf[2]=coord_constraint.at(2).at(1568+i);
		q[count+4+i]=qf;
	}
    
	for (int i=0; i<count+fixed;i++){
	std::cout << "x:" << q[i][0]<< ", y:" << q[i][1] << ", z:" << q[i][2] << std::endl;
	}

    DeformationGraph graph(nb_nodes, nodes, nb_vertices, human, k_nearest);
    graph.print();
    gaussNewton(graph, count, human, q);
    //~ graph.print();

    double* deformed[nb_vertices];
    for (int i = 0 ; i < nb_vertices ; i++)
    {
        double* node = new double[3];
        deformed[i] = node;
    }

    graph.predictAll(nb_vertices, human, deformed);
    //~ std::cout << "deformed:\n";
    //~ for (int i = 0 ; i < nb_vertices ; i++)
        //~ std::cout << deformed[i][0] << " " << deformed[i][1] << " " << deformed[i][2] << "\n";
     
    std::ofstream writer;
	writer.open("export.xyz");
	if (writer.is_open()) 
	{
		for (int i = 0 ; i < nb_vertices ; i++)
			writer << deformed[i][0] << " " << deformed[i][1] << " " << deformed[i][2] << "\n";
	}
	else 
	{
		std::cout<<"could not open writer"<<std::endl;
	 }
	writer.close();
	
}
