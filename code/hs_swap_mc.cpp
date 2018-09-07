#include <map>
#include <set>
#include <list>
#include <cmath>
#include <ctime>
#include <deque>
#include <queue>
#include <stack>
#include <string>
#include <bitset>
#include <cstdio>
#include <limits>
#include <vector>
#include <climits>
#include <cstring>
#include <cstdlib>
#include <fstream>
#include <numeric>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <unordered_map>
#include <tuple>
#include <array>
#include <random>
#include <regex>
#include <set>
#include <omp.h>
#include <chrono>
//#include <windows.h>

//volatile DWORD dwStart;

using namespace std;

const double pi = 3.1415926535897;

//convenient object to describe a particle
struct Particle {
	double rx, ry, rz;
	int type;
	int index;
};

//convenient object to store the details of the simulation
struct State {
	//explicitly loaded in
	int N; //number of particles
	double L; //box edge length
	vector< double > type_to_diameter;
	vector< string > type_to_atom;
	vector< Particle > particles; //stores the particles
	//created via initialization function
	int num_types;
	double V_particles;
	double V;
	double eta;
	double mean_diameter, min_diameter, max_diameter;
	vector< int > type_to_N_type; //maps particle type to the number of them
	//vector< tuple< int > > swap_pair_options; //array of possible particle pairs to swap

	//does any initialization needed after loading in relevant details
	void Initialize(){
		//calculate the volume
		V = pow(L, 3);

		//count the number of each type present
		Particle particle;
		num_types = type_to_diameter.size();
		type_to_N_type.clear();
		type_to_N_type.resize(num_types, 0);
		eta = 0.0;
		V_particles = 0.0;
		mean_diameter = 0.0;
		min_diameter = 1.0e100;
		max_diameter = 0.0;
		for (auto it = particles.begin(); it < particles.end(); it++){
			particle = *it;
			V_particles = V_particles + pi*pow(type_to_diameter[particle.type], 3) / 6.0;
			mean_diameter = mean_diameter + type_to_diameter[particle.type];
			min_diameter = min(min_diameter, type_to_diameter[particle.type]);
			max_diameter = max(max_diameter, type_to_diameter[particle.type]);
			type_to_N_type[particle.type]++;
		}
		mean_diameter = mean_diameter / (double)N;

		//calculate the volume fraction
		eta = V_particles / V;
	}
};

//this is a cell list structure
struct CellList{
	int N_cells = -1;
	double L_cell;
	bool active = false;
	vector< vector< vector< list<int> > > > cells;
	vector< tuple<int, int, int> > neighbor_cell_directions = { { make_tuple(1, 0, 0), make_tuple(0, 1, 0), make_tuple(0, 0, 1),
																  make_tuple(1, 1, 0), make_tuple(1, 0, 1), make_tuple(0, 1, 1),
																  make_tuple(1, 1, 1),
																  make_tuple(1, -1, 0), make_tuple(1, 0, -1), make_tuple(0, 1, -1),
																  make_tuple(-1, 1, 1), make_tuple(1, -1, 1), make_tuple(1, 1, -1)
															  } };

	//structure the cell list
	void PrepareCellList(State &state, bool force_rebuild=false){
		int N_cells_prev = N_cells;
		bool active_prev = active;
		N_cells = (int)(state.L / state.type_to_diameter.back());
		active = N_cells > 3;
		if (active){
			//if not the same number of cells or if they have not been maintained do not rebuild
			if (N_cells != N_cells_prev || !active_prev || force_rebuild){
				cout << "rebuilding cells" << endl;
				//build the cells
				int cell_x, cell_y, cell_z;
				cells.clear();
				L_cell = state.L / (float)N_cells;
				cells.resize(N_cells);
				for (cell_x = 0; cell_x < N_cells; cell_x++)
				{
					cells[cell_x].resize(N_cells);
					for (cell_y = 0; cell_y < N_cells; cell_y++)
					{
						cells[cell_x][cell_y].resize(N_cells);
					}
				}
				//load particles into the cells
				for (auto it = state.particles.begin(); it < state.particles.end(); it++){
					tie(cell_x, cell_y, cell_z) = FindCell((*it));
					cells[cell_x][cell_y][cell_z].push_back(distance(state.particles.begin(), it));
				}
			}
		}
	}

	//finds the cell a particle belongs to (never gets used outside of struct so no "active" check)
	tuple<int, int, int> FindCell(Particle &particle){
		int cell_x = (int)(particle.rx / L_cell);
		int cell_y = (int)(particle.ry / L_cell);
		int cell_z = (int)(particle.rz / L_cell);
		return make_tuple(cell_x, cell_y, cell_z);
	}

	//swaps a particle from one cell to another
	void SingleParticleUpdate(Particle &particle_new, Particle &particle_current){
		if (active){
			int cell_x_new, cell_y_new, cell_z_new;
			int cell_x_current, cell_y_current, cell_z_current;
			tie(cell_x_current, cell_y_current, cell_z_current) = FindCell(particle_current);
			//find the particle in it
			auto it = find(cells[cell_x_current][cell_y_current][cell_z_current].begin(),
				cells[cell_x_current][cell_y_current][cell_z_current].end(), particle_current.index);
			//check to make sure it is there///////////////////////////////////
			if (it == cells[cell_x_current][cell_y_current][cell_z_current].end()){
				cout << "Failed to find particle " << particle_current.index << " in cell list" << endl; getchar();
			}
			////////////////////////////////////////////////////////////////////
			if (particle_current.index != particle_new.index){
				cout << "Particle mismatch upon swapping cells" << endl; getchar();
			}
			//perform the swap
			cells[cell_x_current][cell_y_current][cell_z_current].erase(it);
			tie(cell_x_new, cell_y_new, cell_z_new) = FindCell(particle_new);
			cells[cell_x_new][cell_y_new][cell_z_new].push_back(particle_new.index);
		}
	}
};

//structure for the volume fractions to be sampled
struct ExtendedState{
	//double eta;
	//double dr_max;
	//double swap_threshold;
	long long int num_trans=1;
	long long int num_trans_accept=0;
	long long int num_swaps=1;
	long long int num_swaps_accept=0;
};

//function definitions
void LoadSimulationData(
	double &eta, long long int &steps_equil, double &hours_equil,
	long long int &steps_prod, double &hours_prod,
	int &skip_steps, int &skip_steps_traj, long long int &skip_steps_check, double &acc_target_translation, double &acc_target_swap, double &frac_trans, int &seed);
//void LoadExtendedStateData(ExtendedState &extended_state);
bool CheckParticleOverlapRSA(State &state, Particle &particle);
bool NearestNeighborOverlap(State &state, Particle &particle1, Particle &particle2);
bool CheckParticleOverlap(State &state, CellList &cell_list, Particle &particle);
bool CheckAllParticleOverlaps(State &state, CellList &cell_list);
void Compress(State &state, CellList &cell_list, double eta_tgt, double scale, double dr_max, double swap_threshold, int seed);
void RandomSequentialAddition(State &state, int seed, int max_attempts = 100000);
void WriteConfig(string file, State &state);
void WriteState(string file, State &state);
bool ReadState(string file, State &state);
tuple<bool, bool> AttemptParticleSwap(State &state, CellList &cell_list, int index_1, int index_2, double swap_threshold);
bool AttemptParticleTranslation(State &state, CellList &cell_list, int index, double dx, double dy, double dz);
void WriteAcceptanceStats(ofstream &acceptance_stats_output, ExtendedState &extended_state);
void MonteCarlo(State &state, CellList &cell_list, ExtendedState &extended_state, long long int steps, double hours, int skip_steps, int skip_steps_traj, long long int skip_steps_check,
	string simulation_name, double dr_max, double swap_threshold, double frac_trans, int seed);
double TranslationTuner(double acc_target, int num_steps_for_stats, int num_trials,
												State &state, CellList &cell_list, ExtendedState &extended_state, double swap_threshold, double frac_trans);
double SwapTuner(double acc_target, int num_steps_for_stats, int num_trials,
								 State &state, CellList &cell_list, ExtendedState &extended_state, double dr_max, double frac_trans);

//this will drive everything based on command line input values
int main(){
	//the simulation state
	State state;
	State state_new;
	CellList cell_list;
	CellList cell_list_new;
	ExtendedState extended_state;

	//load in the relevant simulation details
	double eta;
	long long int steps_equil, steps_prod, skip_steps_check;
	double hours_equil, hours_prod;
	int skip_steps, skip_steps_traj, seed;
	double acc_target_translation, acc_target_swap, frac_trans;
	bool read_state;
	LoadSimulationData(eta, steps_equil, hours_equil, steps_prod,
		hours_prod, skip_steps, skip_steps_traj, skip_steps_check, acc_target_translation, acc_target_swap, frac_trans, seed); cout << endl;

	//read in an initial state or perform random sequential addition
	read_state = ReadState("state_init.txt", state); cout << endl;
	if (!read_state){
		RandomSequentialAddition(state, seed);
		WriteState("state_rsa.txt", state);
	}
	state.Initialize();
	cell_list.PrepareCellList(state, true);

	//load in the extended states data
	//LoadExtendedStateData(extended_state);


	//set the parameters for tuning acceptance rates
	int num_trials = 20;
	double dr_max, swap_threshold;


	//compress the thing if needed to achieve target and set active extended state
	if (state.eta < eta){
		cout << endl;
		cout << "Tuning before compression..." << endl;
		//tune the translation assuming everything swaps
		dr_max = TranslationTuner(acc_target_translation, 1000000, num_trials,
															 state, cell_list, extended_state, swap_threshold=100.0, frac_trans=1.0);
		cout << endl;
	  //tune the swap with just a few translation moves
		swap_threshold = SwapTuner(acc_target_swap, 3000000, num_trials,
															 state, cell_list, extended_state, dr_max=dr_max, frac_trans=0.40);
		Compress(state, cell_list, eta, 0.999, dr_max, swap_threshold, seed);
	}
	state.Initialize();
	cell_list.PrepareCellList(state, true);



	//display number of cells and mean particle diameter
	cout << endl;
	cout << "Number of cells: " << cell_list.N_cells << endl;
	cout << "Mean particle diameter: " << state.mean_diameter << endl;
	cout << "Min particle diameter: " << state.min_diameter << endl;
	cout << "Max particle diameter: " << state.max_diameter << endl << endl;
	cout << endl;

	/////////////////////////////////////////////////////////////
	//optimize the tunable sampling variables before actual runs
	cout << endl;
	cout << "Tuning before actual simulation..." << endl;
	//tune the translation assuming everything swaps
	dr_max = TranslationTuner(acc_target_translation, 1000000, num_trials,
														 state, cell_list, extended_state, swap_threshold=100.0, frac_trans=1.0);
	cout << endl;
  //tune the swap with just a few translation moves
	swap_threshold = SwapTuner(acc_target_swap, 3000000, num_trials,
														 state, cell_list, extended_state, dr_max=dr_max, frac_trans=0.40);
	cout << endl;
	/////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////


	//perform monte carlo iterations
	cout << endl;
	WriteState("state_pre_equil.txt", state);
	MonteCarlo(state, cell_list, extended_state, steps_equil, hours_equil, skip_steps, skip_steps_traj,
		skip_steps_check, "equilibration", dr_max, swap_threshold, frac_trans, seed+1);

	cout << endl;
	WriteState("state_pre_prod.txt", state);
	MonteCarlo(state, cell_list, extended_state, steps_prod, hours_prod, skip_steps, skip_steps_traj,
		skip_steps_check, "production", dr_max, swap_threshold, frac_trans, seed);
	WriteState("state_post_prod.txt", state);

	//write one final config
	WriteConfig("trajectory_final.xyz", state);

	return 0;
}

//reads in the intial details like volume fraction and particle composition
void LoadSimulationData(
	double &eta, long long int &steps_equil, double &hours_equil,
	long long int &steps_prod, double &hours_prod,
	int &skip_steps, int &skip_steps_traj, long long int &skip_steps_check, double &acc_target_translation, double &acc_target_swap, double &frac_trans, int &seed){

	ifstream in_file;
	string line;
	in_file.open("./simulation.txt");
	if (!in_file) { cerr << "Unable to open simulation.txt"; exit(1); }

	//prepare the regex's and match
	regex num_re("([0-9\\.\\-]+)");
	smatch match;

	//read in the eta
	getline(in_file, line);
	regex_search(line, match, num_re);
	eta = stod(match.str(1));

	//read in equil details
	getline(in_file, line);
	regex_search(line, match, num_re);
	steps_equil = stoll(match.str(1));
	getline(in_file, line);
	regex_search(line, match, num_re);
	hours_equil = stod(match.str(1));

	//read in prod details
	getline(in_file, line);
	regex_search(line, match, num_re);
	steps_prod = stoll(match.str(1));
	getline(in_file, line);
	regex_search(line, match, num_re);
	hours_prod = stod(match.str(1));

	//read in skip_steps
	getline(in_file, line);
	regex_search(line, match, num_re);
	skip_steps = stoi(match.str(1));

	//read in skip_steps_traj
	getline(in_file, line);
	regex_search(line, match, num_re);
	skip_steps_traj = stoi(match.str(1));

	//read in skip_steps_check
	getline(in_file, line);
	regex_search(line, match, num_re);
	skip_steps_check = stoll(match.str(1));

	//read in acc_target_translation
	getline(in_file, line);
	regex_search(line, match, num_re);
	acc_target_translation = stod(match.str(1));

	//read in acc_target_swap
	getline(in_file, line);
	regex_search(line, match, num_re);
	acc_target_swap = stod(match.str(1));

	//read in frac_trans
	getline(in_file, line);
	regex_search(line, match, num_re);
	frac_trans = stod(match.str(1));

	//read in seed
	getline(in_file, line);
	regex_search(line, match, num_re);
	seed = stoi(match.str(1));

	//close the file
	in_file.close();

	cout << "Simulation details:" << endl;
	cout << "eta = " << eta << endl;
	cout << "steps_equil = " << steps_equil << endl;
	cout << "hours_equil = " << hours_equil << endl;
	cout << "steps_prod = " << steps_prod << endl;
	cout << "hours_prod = " << hours_prod << endl;
	cout << "skip_steps = " << skip_steps << endl;
	cout << "skip_steps_traj = " << skip_steps_traj << endl;
	cout << "skip_steps_check = " << skip_steps_check << endl;
	cout << "acc_target_translation = " << acc_target_translation << endl;
	cout << "acc_target_swap = " << acc_target_swap << endl;
	cout << "frac_trans = " << frac_trans << endl;
	cout << "seed = " << seed << endl;
	cout << endl;
}

//reads in intialization data for the run
/*void LoadExtendedStateData(ExtendedState &extended_state){
	ifstream in_file;
	string line;
	in_file.open("./extended_state.txt");
	if (!in_file) { cerr << "Unable to open extended_state.txt"; exit(1); }

	//prepare the regex's and match
	regex num_re("([0-9\\.]+)");
	sregex_iterator next, end;
	smatch match;

	//read in eta
	getline(in_file, line);
	regex_search(line, match, num_re);
	extended_state.eta = stod(match.str(1));
	extended_state.num_trans_accept = 0;
	extended_state.num_trans = 0;
	extended_state.num_swaps_accept = 0;
	extended_state.num_swaps= 0;

	//read in dr_max
	getline(in_file, line);
	regex_search(line, match, num_re);
	extended_state.dr_max = stod(match.str(1));

	//read in swap_threshold
	getline(in_file, line);
	regex_search(line, match, num_re);
	extended_state.swap_threshold = stod(match.str(1));

	cout << "Loaded extended state data" << endl;
	cout << "eta = " << extended_state.eta << endl;
	cout << "dr_max = " << extended_state.dr_max << endl;
	cout << "swap_threshold = " << extended_state.swap_threshold << endl << endl;

	in_file.close();
}*/

//checks if a specified particle overlaps with another via brute force methods for RSA only!
bool CheckParticleOverlapRSA(State &state, Particle &particle){
	Particle existing_particle;
	double squared_distance;
	double min_squared_distance;
	double dx, dy, dz;

	//iterate over existing particles to see if overlap
	for (auto it = state.particles.begin(); it < state.particles.end(); it++){
		existing_particle = *it;
		dx = particle.rx - existing_particle.rx;
		dy = particle.ry - existing_particle.ry;
		dz = particle.rz - existing_particle.rz;
		dx = dx - round(dx / state.L) * state.L;
		dy = dy - round(dy / state.L) * state.L;
		dz = dz - round(dz / state.L) * state.L;

		squared_distance = pow(dx, 2) + pow(dy, 2) + pow(dz, 2);

		//check for overlap and if found push iterator to end
		min_squared_distance = pow(state.type_to_diameter[particle.type] + state.type_to_diameter[existing_particle.type], 2) / 4.0;
		if (squared_distance < min_squared_distance){
			if (particle.index != existing_particle.index){
				return false;
			}
		}
	}
	return true;
}

//checks for a nearest neighbor overlap between two particles and will work for the edge case of the same particle
bool NearestNeighborOverlap(State &state, Particle &particle1, Particle &particle2){
	if (particle1.index == particle2.index)
		return false;
	double dx, dy, dz, squared_distance, min_squared_distance;
	dx = particle1.rx - particle2.rx;
	dy = particle1.ry - particle2.ry;
	dz = particle1.rz - particle2.rz;
	dx = dx - round(dx / state.L) * state.L;
	dy = dy - round(dy / state.L) * state.L;
	dz = dz - round(dz / state.L) * state.L;
	squared_distance = pow(dx, 2) + pow(dy, 2) + pow(dz, 2);
	min_squared_distance = pow(state.type_to_diameter[particle1.type] + state.type_to_diameter[particle2.type], 2) / 4.0;
	if (squared_distance < min_squared_distance)
		return true;
	else
		return false;
}

//checks if a specified particle overlaps with another via brute force methods
bool CheckParticleOverlap(State &state, CellList &cell_list, Particle &particle){
	Particle existing_particle;

	//check the cell list to see if active, if not just brute force it
	if (cell_list.active){
		//get cell that the particle resides in
		int cell_x_part, cell_y_part, cell_z_part;
		tie(cell_x_part, cell_y_part, cell_z_part) = cell_list.FindCell(particle);

		//loop over neighboring cells
		//int cell_x, cell_y, cell_z;
		int cell_x_wrpd, cell_y_wrpd, cell_z_wrpd;
		for (int cell_x = cell_x_part - 1; cell_x <= cell_x_part + 1; cell_x++){
			for (int cell_y = cell_y_part - 1; cell_y <= cell_y_part + 1; cell_y++){
				for (int cell_z = cell_z_part - 1; cell_z <= cell_z_part + 1; cell_z++){
					cell_x_wrpd = cell_x + (int)(cell_x == -1)*cell_list.N_cells - (int)(cell_x == cell_list.N_cells)*cell_list.N_cells;
					cell_y_wrpd = cell_y + (int)(cell_y == -1)*cell_list.N_cells - (int)(cell_y == cell_list.N_cells)*cell_list.N_cells;
					cell_z_wrpd = cell_z + (int)(cell_z == -1)*cell_list.N_cells - (int)(cell_z == cell_list.N_cells)*cell_list.N_cells;
					//loop over particles
					for (auto it = cell_list.cells[cell_x_wrpd][cell_y_wrpd][cell_z_wrpd].begin(); it != cell_list.cells[cell_x_wrpd][cell_y_wrpd][cell_z_wrpd].end(); it++){
						existing_particle = state.particles[(*it)];
						if (NearestNeighborOverlap(state, particle, existing_particle))
							return false;
					}
				}
			}
		}
		return true;
	}
	//the brute force option if all else fails
	else{
		for (auto it = state.particles.begin(); it < state.particles.end(); it++){
			existing_particle = *it;
			if (NearestNeighborOverlap(state, particle, existing_particle))
				return false;
		}
		return true;
	}
}

//checks for any possible overlap
bool CheckAllParticleOverlaps(State &state, CellList &cell_list){
	Particle particle1, particle2;
	int d_cell_x, d_cell_y, d_cell_z;
	int cell_x_n, cell_y_n, cell_z_n;

	//check the cell list to see if active, if not just brute force it
	if (cell_list.active){
		for (int cell_x = 0; cell_x < cell_list.N_cells; cell_x++){
			for (int cell_y = 0; cell_y < cell_list.N_cells; cell_y++){
				for (int cell_z = 0; cell_z < cell_list.N_cells; cell_z++){
					//check intra cell first
					for (auto iter1 = cell_list.cells[cell_x][cell_y][cell_z].begin(); iter1 != cell_list.cells[cell_x][cell_y][cell_z].end(); iter1++){
						auto iter2 = iter1; advance(iter2, 1);
						for (; iter2 != cell_list.cells[cell_x][cell_y][cell_z].end(); iter2++){
							particle1 = state.particles[(*iter1)];
							particle2 = state.particles[(*iter2)];
							if (NearestNeighborOverlap(state, particle1, particle2))
								return false;
						}
					}

					//check with neighbor cells now
					for (auto iter = cell_list.neighbor_cell_directions.begin(); iter < cell_list.neighbor_cell_directions.end(); iter++)
					{
						//extract displacements to find 1/2 of total neighbors
						tie(d_cell_x, d_cell_y, d_cell_z) = *iter;
						cell_x_n = cell_x + d_cell_x;
						cell_y_n = cell_y + d_cell_y;
						cell_z_n = cell_z + d_cell_z;

						//wrap the cells using periodic boundary conditions
						cell_x_n = cell_x_n + (int)(cell_x_n == -1)*cell_list.N_cells - (int)(cell_x_n == cell_list.N_cells)*cell_list.N_cells;
						cell_y_n = cell_y_n + (int)(cell_y_n == -1)*cell_list.N_cells - (int)(cell_y_n == cell_list.N_cells)*cell_list.N_cells;
						cell_z_n = cell_z_n + (int)(cell_z_n == -1)*cell_list.N_cells - (int)(cell_z_n == cell_list.N_cells)*cell_list.N_cells;

						//loop over the particles in each cell
						for (auto iter1 = cell_list.cells[cell_x][cell_y][cell_z].begin(); iter1 != cell_list.cells[cell_x][cell_y][cell_z].end(); iter1++){
							for (auto iter2 = cell_list.cells[cell_x_n][cell_y_n][cell_z_n].begin(); iter2 != cell_list.cells[cell_x_n][cell_y_n][cell_z_n].end(); iter2++){
								particle1 = state.particles[(*iter1)];
								particle2 = state.particles[(*iter2)];
								if (NearestNeighborOverlap(state, particle1, particle2))
									return false;
							}
						}
					}
				}
			}
		}
		return true;
	}
	//the brute force N^2 option if all else fails
	else{
		for (auto iter1 = state.particles.begin(); iter1 < state.particles.end(); iter1++){
			for (auto iter2 = iter1 + 1; iter2 < state.particles.end(); iter2++){
				particle1 = *iter1;
				particle2 = *iter2;
				if (NearestNeighborOverlap(state, particle1, particle2))
					return false;
			}
		}
		return true;
	}
}

//compresses the system to a target volume fraction
void Compress(State &state, CellList &cell_list, double eta_tgt, double scale, double dr_max, double swap_threshold, int seed){
	mt19937_64 rng_index, rng_x, rng_y, rng_z;
	mt19937_64 rng_index1, rng_index2;
	rng_index.seed(99991+seed); rng_x.seed(3465+seed); rng_y.seed(1503+seed); rng_z.seed(3387773+seed);
	rng_index1.seed(109 + seed); rng_index2.seed(1814 + seed);
	uniform_int_distribution<int> r_index(0, state.N - 1);
	uniform_real_distribution<double> r_dr(-dr_max, dr_max);
	double V_tgt = state.V_particles / eta_tgt;
	double V_new, L_new;
	double scale_L;
	bool target_reached = false;

	cout << "Starting compression from " << state.eta << " to " << eta_tgt << "..." << endl;

	//loop until the target eta is reached
	while (!target_reached){
		V_new = state.V*scale;
		L_new = pow(V_new, 1.0 / 3.0);
		state.V = V_new;
		scale_L = L_new / state.L;
		state.L = L_new;
		state.eta = state.V_particles / state.V;

		//move particles a bit
		for (auto it = state.particles.begin(); it < state.particles.end(); it++){
			(*it).rx = (*it).rx * scale_L;
			(*it).rx = (*it).rx - floor((*it).rx / state.L)*state.L;
			(*it).ry = (*it).ry * scale_L;
			(*it).ry = (*it).ry - floor((*it).ry / state.L)*state.L;
			(*it).rz = (*it).rz * scale_L;
			(*it).rz = (*it).rz - floor((*it).rz / state.L)*state.L;
		}

		//initialize the state and cell list
		state.Initialize();
		cell_list.PrepareCellList(state, true);

		cout << endl << "Current eta: " << state.eta << endl;
		cout << "Current num_cells: " << cell_list.N_cells << endl << endl;

		//check for overlaps
		while (!CheckAllParticleOverlaps(state, cell_list)){
			//perform a bunch of MC moves to hopefully free up overlaps
			for (int i = 0; i < 100 * state.N; i++){
				AttemptParticleTranslation(state, cell_list, r_index(rng_index), r_dr(rng_x), r_dr(rng_y), r_dr(rng_z));
				AttemptParticleSwap(state, cell_list, r_index(rng_index1), r_index(rng_index2), swap_threshold);
			}
		}

		//check if done and release compression a bit to exact desired density
		if (state.V < V_tgt){
			L_new = pow(V_tgt, 1.0 / 3.0);
			state.V = V_tgt;
			state.L = L_new;
			state.eta = state.V_particles / state.V;

			cout << "Finished compression to a volume fraction of " << state.eta << endl << endl;
			target_reached = true;
		}
		WriteConfig("trajectory_compress.xyz", state);
	}
}

//general rsa based on an input file
void RandomSequentialAddition(State &state, int seed, int max_attempts){
	//temporary storage
	int n;
	string atom;
	double diameter;

	//for reading the file
	ifstream in_file;
	string line;
	in_file.open("./composition.txt");
	if (!in_file) { cerr << "Unable to open composition.txt"; exit(1); }

	//prepare the regex's and match
	regex num_re("([0-9\\.]+)");
	regex particle_info_re("\\(([a-zA-Z]+)\\s*,\\s*([0-9\\.]+)\\s*,\\s*([0-9]+)\\)");
	sregex_iterator next, end;
	smatch match;

	//read in the volume fraction
	getline(in_file, line);
	regex_search(line, match, num_re);
	double eta = stod(match.str(1));

	//read in the particles and clear out state
	getline(in_file, line);
	next = sregex_iterator(line.begin(), line.end(), particle_info_re);
	end = sregex_iterator();
	state.N = 0;
	state.L = 0.0;
	state.V_particles = 0.0;
	state.particles.clear();
	state.type_to_atom.clear();
	state.type_to_diameter.clear();
	state.type_to_N_type.clear();

	cout << "Reading in the composition data" << endl;

	//fetch the composition data
	while (next != end) {
		match = *next;
		atom = match.str(1);
		diameter = stod(match.str(2));
		n = stoi(match.str(3));
		state.type_to_atom.push_back(atom);
		state.type_to_diameter.push_back(diameter);
		state.type_to_N_type.push_back(n);
		state.N = state.N + n;
		state.V_particles = state.V_particles + n*(pi / 6.0)*pow(diameter, 3.0);
		next++;
	}
	in_file.close();

	//calculate the L needed for the desired volume fraction
	state.L = pow(state.V_particles / eta, 1.0 / 3.0);


	//random number generator
	mt19937_64 rng_x, rng_y, rng_z;
	rng_x.seed(12+seed); rng_y.seed(100+seed); rng_z.seed(3+seed);
	uniform_real_distribution<double> unif(0, state.L);

	//items for attempting insertions
	//int attempt, n;
	//bool overlap_free;
	//int particle_type;
	//Particle particle;

	cout << "Starting RSA generation at eta = " << eta << endl;


	//loop over the types and how many of them
	Particle particle;
	bool overlap_free;
	int attempt = 0;
	int index = 0;
	for (int i = 0; i < (int)state.type_to_N_type.size() && attempt < max_attempts; i++){

		//loop over attempts to insert desired number of each type
		attempt = 0; n = 0;
		while (attempt < max_attempts && n < state.type_to_N_type[i]){
			particle.rx = unif(rng_x); particle.ry = unif(rng_y); particle.rz = unif(rng_z);
			particle.type = i;
			particle.index = index;

			//check for an overlap
			overlap_free = CheckParticleOverlapRSA(state, particle);

			//add new particle if possible
			if (overlap_free){
				attempt = 0; n++; state.particles.push_back(particle);
				index++;
			}
			else{
				attempt++;
			}
		}
	}

	//provide some feedback
	if (attempt < max_attempts)
		cout << "RSA completed!" << endl << endl;

}

//take the state and write out an xyz file of it
void WriteConfig(string file, State &state){
	ofstream output(file, ios::app);
	output << state.N << endl;
	output << std::setprecision(10) << " " << "L = " << state.L << endl;

	Particle particle;
	string atom;
	for (auto it = state.particles.begin(); it < state.particles.end(); it++){
		particle = *it;
		atom = state.type_to_atom[particle.type] + "       ";
		output << std::setprecision(10) << "  " << atom << particle.rx << "       " << particle.ry << "       " << particle.rz << "       " << state.type_to_diameter[particle.type] / 2.0 << endl;
	}
	output.close();
}

//serialize the state object
void WriteState(string file, State &state){
	ofstream output(file, ios::app);
	output << state.N << " " << state.L << endl; //save the N and L
	for (auto it = state.type_to_atom.begin(); it < state.type_to_atom.end(); it++){ //save the map from type to atom
		output << "(" << distance(state.type_to_atom.begin(), it) << "," << *it << ")" << " ";
	}
	output << endl;
	for (auto it = state.type_to_diameter.begin(); it < state.type_to_diameter.end(); it++){ //save the map from type to atom
		output << std::setprecision(10) << "(" << distance(state.type_to_diameter.begin(), it) << "," << *it << ")" << " ";
	}
	output << endl;
	Particle particle;
	//string atom_type;
	for (auto it = state.particles.begin(); it < state.particles.end(); it++){
		particle = *it;
		output << std::setprecision(10) << "(" << particle.rx << "," << particle.ry << "," << particle.rz << "," << particle.type << "," << particle.index << ")" << " ";
	}
	output << endl;
	output.close();
}

//read in the serialized state object
bool ReadState(string file, State &state){
	ifstream in_file;
	string line;
	in_file.open(file);
	if (!in_file) {cerr << "Unable to open state file.txt"; return false/*exit(1)*/;}

	//prepare the regex's and match
	regex N_and_L_re("([0-9]+)\\s+([0-9\\.]+)");
	regex type_and_atom_re("\\(([0-9]+),([a-zA-Z]+)\\)");
	regex type_and_diameter_re("\\(([0-9]+),([0-9\\.]+)\\)");
	regex particles_re("\\(([0-9\\.eE\\-]+),([0-9\\.eE\\-]+),([0-9\\.eE\\-]+),([0-9]+),([0-9]+)\\)");
	sregex_iterator next, end;
	smatch match;

	//read in N and L
	getline(in_file, line);
	regex_search(line, match, N_and_L_re);
	state.N = stoi(match.str(1)); state.L = stod(match.str(2));

	//read in the type to atom
	getline(in_file, line);
	//regex_search(line, match, type_and_atom_re);
	next = sregex_iterator(line.begin(), line.end(), type_and_atom_re);
	end = sregex_iterator();
	state.type_to_atom.clear();
	while (next != end) {
		match = *next;
		cout << match.str(1) << "," << match.str(2) << "\n";
		//state.type_to_atom[stoi(match.str(1))] = match.str(2);
		state.type_to_atom.push_back(match.str(2));
		next++;
	}

	//read in the type to diameter
	getline(in_file, line);
	next = sregex_iterator(line.begin(), line.end(), type_and_diameter_re);
	end = sregex_iterator();
	state.type_to_diameter.clear();
	while (next != end) {
		match = *next;
		cout << match.str(1) << "," << match.str(2) << "\n";
		//state.type_to_diameter[stoi(match.str(1))] = stod(match.str(2));
		state.type_to_diameter.push_back(stod(match.str(2)));
		next++;
	}

	//read in particles
	getline(in_file, line);
	next = sregex_iterator(line.begin(), line.end(), particles_re);
	end = sregex_iterator();
	Particle particle;
	state.particles.clear();
	while (next != end) {
		match = *next;
		particle.rx = stod(match.str(1)); particle.ry = stod(match.str(2)); particle.rz = stod(match.str(3));
		particle.type = stoi(match.str(4));
		particle.index = stoi(match.str(5));
		state.particles.push_back(particle);
		next++;
	}

	cout << endl <<  "Read in " << state.particles.size() << " particles" << endl;
	return true;
}

//attempt a particle swap
tuple<bool, bool> AttemptParticleSwap(State &state, CellList &cell_list, int index_1, int index_2, double swap_threshold){
	//save the types
	int type_1_curr = state.particles[index_1].type;
	int type_2_curr = state.particles[index_2].type;

	//if the same type just do it
	if (type_1_curr == type_2_curr)
		return make_tuple(true, false); //return true;

	//only swap similar particles that deviate
	double d1 = state.type_to_diameter[type_1_curr];
	double d2 = state.type_to_diameter[type_2_curr];
	double diff = abs(d1 - d2);
	if (diff > swap_threshold*(state.max_diameter-state.min_diameter))
		return make_tuple(false, false); //return false;

	//change the types if different
	state.particles[index_1].type = type_2_curr;
	state.particles[index_2].type = type_1_curr;

	//check for overlaps (only for the larger particle that is)
	if (d1 > d2){
		if (CheckParticleOverlap(state, cell_list, state.particles[index_2]))
			return make_tuple(true, true); //return true;
	}
	else{
		if (CheckParticleOverlap(state, cell_list, state.particles[index_1]))
			return make_tuple(true, true); //return true;
	}
	state.particles[index_1].type = type_1_curr;
	state.particles[index_2].type = type_2_curr;
	return make_tuple(false, true); //return false;

}

//attempts to translate a particle with periodic wrapping
bool AttemptParticleTranslation(State &state, CellList &cell_list, int index, double dx, double dy, double dz){
	Particle particle_translated = state.particles[index];
	particle_translated.rx = particle_translated.rx + dx;
	particle_translated.rx = particle_translated.rx - floor(particle_translated.rx / state.L)*state.L;
	particle_translated.ry = particle_translated.ry + dy;
	particle_translated.ry = particle_translated.ry - floor(particle_translated.ry / state.L)*state.L;
	particle_translated.rz = particle_translated.rz + dz;
	particle_translated.rz = particle_translated.rz - floor(particle_translated.rz / state.L)*state.L;
	//check if overlap and make move if possible
	if (CheckParticleOverlap(state, cell_list, particle_translated)){
		cell_list.SingleParticleUpdate(particle_translated, state.particles[index]);
		state.particles[particle_translated.index] = particle_translated;
		return true;
	}
	else
		return false;
}

//writes out the acceptance stats for each eta
void WriteAcceptanceStats(ofstream &acceptance_stats_output, ExtendedState &extended_state){
	acceptance_stats_output << "Trans: " << (double)extended_state.num_trans_accept / (double)extended_state.num_trans
	<< "    |    Swap: " << (double)extended_state.num_swaps_accept / (double)extended_state.num_swaps << endl;
}

//perform monte carlo steps
void MonteCarlo(State &state, CellList &cell_list, ExtendedState &extended_state, long long int steps, double hours, int skip_steps, int skip_steps_traj, long long int skip_steps_check,
				string simulation_name, double dr_max, double swap_threshold, double frac_trans, int seed)
{
	//random number generator
	mt19937_64 rng_move_type, rng_index, rng_swap_index_1, rng_swap_index_2, rng_x, rng_y, rng_z, rng_eta_change, rng_eta_change_accept;
	rng_move_type.seed(1 + seed); rng_index.seed(991 + seed); rng_swap_index_1.seed(62 + seed); rng_swap_index_2.seed(65781 + seed); rng_x.seed(34 + seed);
	rng_y.seed(103+seed); rng_z.seed(333+seed); rng_eta_change.seed(6+seed); rng_eta_change_accept.seed(600+seed);
	uniform_int_distribution<int> r_index(0, state.N - 1);
	//uniform_real_distribution<double> r_dr(-dr_max, dr_max);
	uniform_real_distribution<double> r_dr(-1.0, 1.0);
	uniform_int_distribution<int> r_eta_change(0, 1); //0=shrink, 1=grow
	uniform_real_distribution<double> r_accept(0, 1);
	bool translation_status, swap_status, count_swap_attempt;
	double move_type;
	string filename, check_filename;
	int min_left;
	int skip_for_stats = 10000000;
	int random_eta_change, active_current;
	double r_dr_x, r_dr_y, r_dr_z;
	int index_1, index_2;

	//open the file fro writing eta statistics
	cout << "Performing " << simulation_name << "..." << endl;

	//zero out the acceptance statistics
	extended_state.num_trans = 1;
	extended_state.num_trans_accept = 0;
	extended_state.num_swaps = 1;
	extended_state.num_swaps_accept = 0;

	//if using time based kill fetch end time
	int minutes = max(0, (int)(hours * 60.0));
	auto end_time = chrono::system_clock::now() + chrono::minutes(minutes);

	long long int i = 1;
	bool iterate = true;
	while (iterate){
		//choose a random type of move
		move_type = r_accept(rng_move_type);

		//random translation move
		if (move_type <= frac_trans){
			r_dr_x = dr_max * r_dr(rng_x);
			r_dr_y = dr_max * r_dr(rng_y);
			r_dr_z = dr_max * r_dr(rng_z);
			translation_status = AttemptParticleTranslation(state, cell_list, r_index(rng_index), r_dr_x, r_dr_y, r_dr_z);
			extended_state.num_trans = extended_state.num_trans + 1;
			extended_state.num_trans_accept = extended_state.num_trans_accept + (int)translation_status;
		}
		//random particle swap
		else{
			index_1 = r_index(rng_swap_index_1);
			index_2 = r_index(rng_swap_index_2);
			tie(swap_status, count_swap_attempt) = AttemptParticleSwap(state, cell_list, index_1, index_2, swap_threshold);
			if(count_swap_attempt){
				extended_state.num_swaps = extended_state.num_swaps + 1;
				extended_state.num_swaps_accept = extended_state.num_swaps_accept + (int)swap_status;
			}
		}

		//message user
		if (i % skip_for_stats == 0){
			cout << "Completed " << i << " steps" << endl;
		}

		//write out trajectory data
		if (i % skip_steps_traj == 0){
			if (simulation_name == "production"){
				cout << "Writing snapshot. " << endl;
				WriteConfig("trajectory.xyz", state);
			}
		}

		//write out checkpoint (only for production)
		if (i % skip_steps_check == 0){
			if (simulation_name == "production"){
				cout << "Writing checkpoint. " << endl;
				check_filename = "state_prod_check__" + to_string(i) + ".txt";
				WriteState(check_filename, state);
			}
		}

		//write out trajectory stats
		if (i % skip_steps == 0){
			filename = simulation_name + "__type_stats.dat";
			min_left = std::chrono::duration_cast<std::chrono::minutes>(end_time - chrono::system_clock::now()).count();
			ofstream acceptance_stats_output(simulation_name + "__acceptance_stats.dat", ios::trunc);
			WriteAcceptanceStats(acceptance_stats_output, extended_state);
			acceptance_stats_output.close();

		}

		//update step count
		i = i + 1;

		//check time our count number
		if (minutes > 0)
			iterate = chrono::system_clock::now() < end_time;
		else
			iterate = i <= steps;
	}

}

//tunes the translation acceptance
double TranslationTuner(double acc_target, int num_steps_for_stats, int num_trials,
												State &state, CellList &cell_list, ExtendedState &extended_state, double swap_threshold, double frac_trans){

	//tune the translation acceptance
	cout << endl << "Tuning the trans acceptance..." << endl;
	double acc;
	double dr_max_max = 5.0*state.mean_diameter;
	double dr_max_min = 0.0;
	double dr_max = (dr_max_min + dr_max_max)/2.0;

	for(int i=0; i<num_trials; i++){
		MonteCarlo(state, cell_list, extended_state, num_steps_for_stats, -1.0, 2*num_steps_for_stats, 2*num_steps_for_stats,
			2*num_steps_for_stats, "tuning_trans", dr_max, swap_threshold, frac_trans, 1021+i);
		acc = (double)extended_state.num_trans_accept/(double)extended_state.num_trans;
		cout << "dr_max: " << dr_max << "    |    acceptance: " << acc << endl;
		cout << "num_trans_accept: " << extended_state.num_trans_accept << "    |    num_trans: " << extended_state.num_trans << endl;
		if (acc < acc_target)
			dr_max_max = dr_max;
		else
			dr_max_min = dr_max;
		dr_max = (dr_max_min + dr_max_max)/2.0;
	}
	return dr_max;
}

//tunes the swap acceptance
double SwapTuner(double acc_target, int num_steps_for_stats, int num_trials,
								 State &state, CellList &cell_list, ExtendedState &extended_state, double dr_max, double frac_trans){

	cout << "Tuning the swap acceptance..." << endl;
	double acc;
	double swap_threshold_max = 1.0;
	double swap_threshold_min = 0.0;
	double swap_threshold = (swap_threshold_min + swap_threshold_max)/2.0;

	for(int i=0; i<num_trials; i++){
		MonteCarlo(state, cell_list, extended_state, num_steps_for_stats, -1.0, 2*num_steps_for_stats, 2*num_steps_for_stats,
			2*num_steps_for_stats, "tuning_swap", dr_max, swap_threshold, frac_trans, 339+i);
		acc = (double)extended_state.num_swaps_accept/(double)extended_state.num_swaps;
		cout << "swap_threshold: " << swap_threshold << "    |    acceptance: " << acc << endl;
		cout << "num_swaps_accept: " << extended_state.num_swaps_accept << "    |    num_swaps: " << extended_state.num_swaps << endl;
		if (acc < acc_target)
			swap_threshold_max = swap_threshold;
		else
			swap_threshold_min = swap_threshold;
		swap_threshold = (swap_threshold_min + swap_threshold_max)/2.0;
	}
	return swap_threshold;
}
