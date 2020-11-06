#include <Klampt/Interface/WorldViewProgram.h>
#include <stdio.h>
#include <iostream>

int main(int argc,const char** argv)
{
	RobotWorld world;
	WorldViewProgram wvp(&world);
	if(!wvp.LoadCommandLine(argc,argv)) {
		return 1;
	}
	auto robot = world.robots[0];
	wvp.Run();
	return 0;
}

void GLUINavigationProgram::Handle_Idle(){

}
