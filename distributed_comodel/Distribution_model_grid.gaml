/**
* Name: Distribution_model
* Distribution model to distribute a thematic model. 
* Author: Lucas Grosjean
* Tags: Visualisation, HPC, distributed ABM, distribution model
*/

model Distribution

import "Continuous_Move_thematic.gaml" as Thematic

global
{
	int end_cycle <- 300;
	int MPI_RANK;
	int MPI_SIZE;
	int grid_width <- 1;
	int grid_height <- 4;
	int size_OLZ <- 5;
	
	float start_thematic_time;
	float end_thematic_time;
	
	file building_shapefile <- file("includes/building.shp");
	//Shape of the environment
	geometry shape <- envelope(building_shapefile); 
	list<people> init_peoples;
	
	init
 	{	
 		start_thematic_time <- machine_time;
 		write("start_thematic_time " + start_thematic_time);
 		
 		create Thematic.Thematic_experiment;
 		create Communication_Agent_MPI;
 		
 		MPI_RANK <- Communication_Agent_MPI[0].MPI_RANK;
 		MPI_SIZE <- Communication_Agent_MPI[0].MPI_SIZE;
		
 		create Partitionning_agent;
 		
 		ask Thematic.Thematic_experiment[0].simulation
 		{
 			init_peoples <- list(people); // get the initial population of people
 		}
 		
		loop current_people over: init_peoples
		{
 			cell c <- cell(current_people.location);
			if(c.rank != MPI_RANK) // kill people not in my cell
			{
				ask current_people
				{
					do die;	
				}
			}
 		}
 	}
 	
 	reflex run_thematic_model // run a cycle of the thematic model
 	{
 		write("distribution step : --------------------------------------" + cycle);
 		ask Thematic.Thematic_experiment[0].simulation
 		{
 			do _step_;
 		}
 	}
 	
 	/*reflex end_distribution_model_no_more_agent // end the distribution model when there are no more agent in the thematic model
 	{
 		ask Thematic.Thematic_experiment[0]
 		{
 			if(length(people) = 0)
 			{
 				end_thematic_time <- machine_time;
 				write("total execution time : " + ((end_thematic_time - start_thematic_time) / 1000) + "second(s)");
 				write("-----------------no more agent to execute-----------------");
 				ask myself
 				{
 					do die;
 				}
 			}
 		}
 	}*/
 	
 	reflex end_distribution_model_end_cycle when: cycle = end_cycle // end the distribution model when we reach the cycle end_cycle
	{
		end_thematic_time <- machine_time;
		write("total execution time : " + ((end_thematic_time - start_thematic_time) / 1000) + "second(s)");
		write("-----------------end_cycle reached-----------------");
		do die;
	}
}

species Partitionning_agent
{
	reflex getCell
	{	
		list<people> peoples;
		list<people> dead_peoples;
		map<int,list<people>> people_to_send;
		
		ask Thematic.Thematic_experiment[0]
		{
			peoples <- people collect each where not dead(each);
			dead_peoples <- people collect each where dead(each);
		}
		loop current_people over: peoples
		{			
			if(not dead(current_people))
			{	
				write("poeple " + current_people.name + " :: " + current_people.location);
				cell c <- cell(current_people.location);
				write("found cell " + c);
				write("found cell rank " + c.rank);
				if(c.rank != MPI_RANK)
				{
					write("yap " + people_to_send[c.rank]);
					if( people_to_send[c.rank] = nil)
					{
						write("nininini");
						people_to_send[c.rank] <- list<people>(current_people);
					}else
					{
						write("not nininini");
						people_to_send[c.rank] << current_people;
					}
					write("yap2 " + people_to_send[c.rank]);
				}
			}	
		}
		do send_people_not_in_my_cell(people_to_send);
	}
	
	action send_people_not_in_my_cell(map<int,list<people>> people_to_send)
	{
		write("sending people : " + people_to_send);
		
		map<int,list<people>> new_people;
		ask Communication_Agent_MPI
		{			
			 new_people <- all_to_all(people_to_send);
		}
		loop peoples over: people_to_send
		{
			ask peoples
			{
				write("killing people " + name);
				do die;
			}
		}
		write("new_people inside " + new_people);
	}
}


grid cell width: grid_width height: grid_height neighbors: 4
{ 
	int rank <- grid_x + (grid_y * grid_width);
	
	list<geometry> OLZ_list;
	map<geometry, int> neighborhood_shape;
	
	geometry OLZ_combined;
	
	/* INNER OLZ */
	geometry OLZ_top_inner <- shape - rectangle(world.shape.width / grid_width, world.shape.height / grid_height) translated_by {0,(size_OLZ / 2),0};
	geometry OLZ_bottom_inner <- shape - rectangle(world.shape.width / grid_width, world.shape.height / grid_height) translated_by {0,-(size_OLZ / 2),0};
	geometry OLZ_left_inner <- shape - rectangle(world.shape.width / grid_width, world.shape.height / grid_height) translated_by {size_OLZ / 2,0,0};
	geometry OLZ_right_inner <- shape - rectangle(world.shape.width / grid_width, world.shape.height / grid_height) translated_by {-(size_OLZ / 2),0,0};
	
	/* CORNER */
	geometry OLZ_bottom_left_inner <- OLZ_left_inner inter OLZ_bottom_inner;
	geometry OLZ_bottom_right_inner <- OLZ_right_inner inter OLZ_bottom_inner;
	geometry OLZ_top_left_inner <- OLZ_left_inner inter OLZ_top_inner;
	geometry OLZ_top_right_inner <- OLZ_right_inner inter OLZ_top_inner;
	
	/* OUTER OLZ */
	geometry OLZ_top_outer <- (shape - rectangle(world.shape.width / grid_width, world.shape.height / grid_height) translated_by {0,(size_OLZ / 2),0}) translated_by {0,-(size_OLZ / 2),0};
	geometry OLZ_bottom_outer <- (shape - rectangle(world.shape.width / grid_width, world.shape.height / grid_height) translated_by {0,-(size_OLZ / 2),0}) translated_by {0,(size_OLZ / 2),0};
	geometry OLZ_left_outer <- (shape - rectangle(world.shape.width / grid_width, world.shape.height / grid_height) translated_by {size_OLZ / 2,0,0}) translated_by {-(size_OLZ / 2),0,0};
	geometry OLZ_right_outer <- (shape - rectangle(world.shape.width / grid_width, world.shape.height / grid_height) translated_by {-(size_OLZ / 2),0,0}) translated_by {(size_OLZ / 2),0,0};
	
	/* ALL INNER OLZ */
	geometry inner_OLZ <- OLZ_top_inner + OLZ_bottom_inner + OLZ_left_inner + OLZ_right_inner;
	
	/* ALL OUTER OLZ */
	geometry outer_OLZ <- OLZ_top_outer + OLZ_bottom_outer + OLZ_left_outer + OLZ_right_outer;
	
	init
	{
		/*write("rank : " + rank);

		write("world.shape.height " + world.shape.height);
		write("world.shape.width " + world.shape.width);*/
		
		// INNER OLZ
		if(grid_y - 1 >= 0)
		{		
			write(""+grid_x + "," + (grid_y-1));
			neighborhood_shape << OLZ_top_inner :: (grid_x + ((grid_y - 1) * grid_width));
			OLZ_combined <- OLZ_combined + OLZ_top_inner;
			OLZ_list << OLZ_top_inner;
		}
		if(grid_y + 1 < grid_height)
		{		
			neighborhood_shape << OLZ_bottom_inner :: (grid_x + ((grid_y + 1) * grid_width));
			OLZ_combined <- OLZ_combined + OLZ_bottom_inner;
			OLZ_list << OLZ_bottom_inner;
		}
		if(grid_x - 1 >=0)
		{		
			neighborhood_shape << OLZ_left_inner :: ((grid_x - 1)  + (grid_y * grid_width));
			OLZ_combined <- OLZ_combined + OLZ_left_inner;
			OLZ_list << OLZ_left_inner;
		}	
		if(grid_x + 1 < grid_width)
		{		
			neighborhood_shape << OLZ_right_inner :: ((grid_x + 1)  + (grid_y * grid_width));
			OLZ_combined <- OLZ_combined + OLZ_right_inner;
			OLZ_list << OLZ_right_inner;
		}
		
		// CORNER
		if(grid_x + 1 < grid_width and grid_y - 1 >= 0)
		{		
			neighborhood_shape << OLZ_top_right_inner :: ((grid_x + 1)  + ((grid_y - 1)  * grid_width));
			OLZ_combined <- OLZ_combined + OLZ_top_right_inner;
			OLZ_list << OLZ_top_right_inner;
		} 
		if(grid_x - 1 >= 0 and grid_y + 1 < grid_height)
		{		
			neighborhood_shape << OLZ_bottom_left_inner :: ((grid_x - 1)  + ((grid_y + 1)  * grid_width));
			OLZ_combined <- OLZ_combined + OLZ_bottom_left_inner;
			OLZ_list << OLZ_bottom_left_inner;
		}
		if(grid_x + 1 < grid_width and grid_y + 1 < grid_height)
		{		
			neighborhood_shape << OLZ_bottom_right_inner :: ((grid_x + 1)  + ((grid_y + 1)  * grid_width));
			OLZ_combined <- OLZ_combined + OLZ_bottom_right_inner;
			OLZ_list << OLZ_bottom_right_inner;
		}
		if(grid_x - 1 >= 0 and grid_y - 1 >= 0)
		{		
			neighborhood_shape << OLZ_top_left_inner :: ((grid_x - 1)  + ((grid_y - 1)  * grid_width));
			OLZ_combined <- OLZ_combined + OLZ_top_left_inner;
			OLZ_list << OLZ_top_left_inner;
		}
	}
	
	aspect default
	{
		draw self.shape color: rgb(#white,125) border:#black;	
		draw "[" + self.grid_x + "," + self.grid_y +"] : " + rank color: #red font:font('Default', 15, #bold) at: {self.location.x, self.location.y};
		
		//if(cell[MPI_RANK] = self)
		//{	
			loop shape_to_display over: neighborhood_shape.keys
			{
				draw shape_to_display color: rgb(#green, 125) border: #black;
			}
			//draw outer_OLZ color: rgb(#red, 125) border: #black;
		//}
	}
}

species Communication_Agent_MPI skills:[MPI_SKILL]
{
 	map<int, unknown> all_to_all(map<int, unknown> data_send)
 	{
 		write("all_to_all from (" + MPI_RANK + ") : " + data_send);
	    map<int, unknown> data_recv <- MPI_ALLTOALL(data_send);
	    write("DATA RECEIVED : " + data_recv);
	    return data_recv;
 	}
}

experiment distribution_experiment type: MPI_EXP  until: (cycle = end_cycle)
{
	list<people> peoples; // trick to print the peoples on the distribution model
	list<building> buildings; // trick to print the buildings on the distribution model
	
	reflex update_peoples
	{
		list<people> peoples_;
		list<building> buildings_;
		ask Thematic.Thematic_experiment[0].simulation 
		{
			 peoples_ <- list(people);
			 buildings_ <- list(building);
		}
		peoples <- peoples_;
		buildings <- buildings_;
	}
	reflex snap
	{
		write("SNAPPING___________________________________ " + cycle);
		int mpi_id <- MPI_RANK;
		//save distribution simulation snapshot
		
		ask simulation
		{	
			save (snapshot("agent")) to: "../output.log/snapshot/" + mpi_id + "/cycle"+ cycle + ".png" rewrite: true;	
		}
	}
	output
	{
		display agent
		{
			agents people value: peoples;	
			agents building value: buildings;
			species cell aspect: default;
		}
	}
}
