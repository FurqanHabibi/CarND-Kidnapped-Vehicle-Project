/*
 * particle_filter.cpp
 *
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang
 */

#include <random>
#include <algorithm>
#include <iostream>
#include <numeric>
#include <math.h> 
#include <iostream>
#include <sstream>
#include <string>
#include <iterator>
#include <limits>
#include <map>

#include "particle_filter.h"
#define EPS 0.0000001 // A very small number

using namespace std;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// cout << "Init" << endl;
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).
	
	// number of particles
	num_particles = 50;

	// setup the gaussian random generator
	default_random_engine gen;
	normal_distribution<double> dist_x(x, std[0]);
	normal_distribution<double> dist_y(y, std[1]);
	normal_distribution<double> dist_theta(theta, std[3]);

	// initialize all particles
	for (int i = 0; i < num_particles; ++i) {
		Particle particle;
		particle.id = i;
		particle.x = dist_x(gen);
		particle.y = dist_y(gen);
		particle.theta = dist_theta(gen);
		particle.weight = 1;
		particles.push_back(particle);
	}

	is_initialized = true;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/
	// cout << "Prediction" << endl;

	default_random_engine gen;

	if ((yaw_rate > -EPS) && (yaw_rate < EPS)) {
		yaw_rate = EPS;
	}
	


	for (int i = 0; i < num_particles; ++i) {
		// predict vehicle's next position and heading based on the velocity and yaw rate
		double x = particles[i].x + 
			(velocity / yaw_rate) * (sin(particles[i].theta + (yaw_rate * delta_t)) - sin(particles[i].theta));
		double y = particles[i].y + 
			(velocity / yaw_rate) * (cos(particles[i].theta) - cos(particles[i].theta + (yaw_rate * delta_t)));
		double theta = particles[i].theta + yaw_rate * delta_t;

		// setup the gaussian random generator
		normal_distribution<double> dist_x(x, std_pos[0]);
		normal_distribution<double> dist_y(y, std_pos[1]);
		normal_distribution<double> dist_theta(theta, std_pos[3]);

		// add random gaussian noise
		particles[i].x = dist_x(gen);
		particles[i].y = dist_y(gen);
		particles[i].theta = dist_theta(gen);
	}
}

void ParticleFilter::dataAssociation(vector<LandmarkObs> predicted, vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.


}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
		const vector<LandmarkObs> &observations, const Map &map_landmarks) {
	// TODO: Update the weights of each particle using a mult-variate Gaussian distribution. You can read
	//   more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
	// NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
	//   according to the MAP'S coordinate system. You will need to transform between the two systems.
	//   Keep in mind that this transformation requires both rotation AND translation (but no scaling).
	//   The following is a good resource for the theory:
	//   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
	//   and the following is a good resource for the actual equation to implement (look at equation 
	//   3.33
	//   http://planning.cs.uiuc.edu/node99.html
	// cout << "updateWeights" << endl;

	// create map structure to cache map_landmarks
	map<int, Map::single_landmark_s> landmarks_map;
	for (const auto& landmarks : map_landmarks.landmark_list) {
		landmarks_map[landmarks.id_i] = landmarks;
	}

	for (int i = 0; i < num_particles; ++i) {
		// cout << "updateWeights - loop " << i << endl;

		// transform the observations to map frame
		vector<LandmarkObs> observations_map;
		for (const auto& obs : observations) {
			LandmarkObs landmarkObs;
			landmarkObs.id = obs.id;
			landmarkObs.x = particles[i].x + (cos(particles[i].theta) * obs.x) - (sin(particles[i].theta) * obs.y);
			landmarkObs.y = particles[i].y + (sin(particles[i].theta) * obs.x) + (cos(particles[i].theta) * obs.y);
			observations_map.push_back(landmarkObs);
		}
		// cout << "updateWeights - loop " << i << " - transformation to map" << endl;

		// associate the observations to map landmarks
		// particles[i].associations.clear();
		// particles[i].sense_x.clear();
		// particles[i].sense_y.clear();
		for (auto& obs : observations_map) {
			double min_distance = numeric_limits<double>::infinity();
			for (const auto& landmark : map_landmarks.landmark_list) {
				double distance = hypot((obs.x - landmark.x_f), (obs.y - landmark.y_f));
				if (distance < min_distance) {
					obs.id = landmark.id_i;
					min_distance = distance;
				}
			}
			// particles[i].associations.push_back(obs.id);
			// particles[i].sense_x.push_back(obs.x);
			// particles[i].sense_y.push_back(obs.y);
		}
		// cout << "updateWeights - loop " << i << " - association" << endl;

		// calculate multivariate gaussian distribution and update weight
		double weight = 1;
		for (const auto& obs : observations_map) {
			// cout << obs.id << obs.x << obs.y << endl;
			weight *= (1 / (2 * M_PI * std_landmark[0] * std_landmark[1])) * 
				exp(-(pow(obs.x - landmarks_map[obs.id].x_f, 2) / (2 * pow(std_landmark[0], 2)))-(pow(obs.y - landmarks_map[obs.id].y_f, 2) / (2 * pow(std_landmark[1], 2))));
		}
		particles[i].weight = weight;
		// cout << "updateWeights - loop " << i << " - updating weight" << endl;
	}

}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
	// cout << "resample" << endl;

	// get the particles' weights in a vector
	weights.clear();
	for (int i = 0; i < num_particles; ++i) {
		weights.push_back(particles[i].weight);
	}

	// setup the discrete distribution generator
	default_random_engine gen;
	discrete_distribution<int> dist(weights.begin(),weights.end());

	// do resampling
	vector<Particle> temp_particles;
	for (int i = 0; i < num_particles; ++i) {
		temp_particles.push_back(particles[dist(gen)]);
	}
	particles = temp_particles;
}

Particle ParticleFilter::SetAssociations(Particle& particle, const vector<int>& associations, 
                                     const vector<double>& sense_x, const vector<double>& sense_y)
{
    //particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
    // associations: The landmark id that goes along with each listed association
    // sense_x: the associations x mapping already converted to world coordinates
    // sense_y: the associations y mapping already converted to world coordinates

    particle.associations= associations;
    particle.sense_x = sense_x;
    particle.sense_y = sense_y;
}

string ParticleFilter::getAssociations(Particle best)
{
	vector<int> v = best.associations;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<int>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseX(Particle best)
{
	vector<double> v = best.sense_x;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseY(Particle best)
{
	vector<double> v = best.sense_y;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
