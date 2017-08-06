/*This is a skeleton code file for use with the Finite Element Method for Problems in Physics.
It uses the deal.II FEM library, dealii.org*/

//Include files
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <unistd.h>
#include <time.h>
#include "FEM_parallel.cc"


using namespace dealii;
void TimeOut(time_t StartTime);

int main (int argc, char *argv[]){

	time_t StartTime;
	StartTime=time(NULL);

	try{
		deallog.depth_console (0);

        loadfiles();

                using namespace dealii;

                Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

                LaplaceProblem<dimension> problemObject(Euler_scheme);
		problemObject.generate_mesh();
		problemObject.setup_system();
                problemObject.solve();
	

	}
	catch (std::exception &exc){
		std::cerr << std::endl << std::endl
		<< "----------------------------------------------------"
		<< std::endl;
		std::cerr << "Exception on processing: " << std::endl
		<< exc.what() << std::endl
		<< "Aborting!" << std::endl
		<< "----------------------------------------------------"
		<< std::endl;

		return 1;
	}
	catch (...){
		std::cerr << std::endl << std::endl
		<< "----------------------------------------------------"
		<< std::endl;
		std::cerr << "Unknown exception" << std::endl
		<< "Aborting!" << std::endl
		<< "----------------------------------------------------"
		<< std::endl;
		return 1;
	}
	TimeOut(StartTime);
	return 0;

}

void TimeOut(time_t StartTime )
{
	int hour, min, sec;
	double duration;
	time_t CurrentTime;

	CurrentTime=time(NULL);
	duration=difftime(CurrentTime, StartTime);
	hour=(int)( duration/3600. );
	min=(int)( (duration - (double)hour*3600.)/60. );
	sec=(int)( duration - (double)hour*3600. - (double)min*60. );
	printf("Time elapsed is %dh %dm %ds \n", hour, min, sec);
}
