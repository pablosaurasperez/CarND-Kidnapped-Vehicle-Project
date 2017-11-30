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

#include "particle_filter.h"

using namespace std;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).
	
	num_particles = 100;
	
	default_random_engine gen;
	double std_x = std[0];
	double std_y = std[1];
	double std_theta = std[2];
	
	//Create Gaussian particles based on first position
	//Gaussian for x
	//normal_distribution<double> dist_x(x, std_x);
	
	//Gaussian for y
	//normal_distribution<double> dist_y(y, std_y);
	
	//Gaussian for theta
	//normal_distribution<double> dist_theta(theta, std_theta);
	
	//Gaussian for x
	normal_distribution<double> dist_x(x, std_x);
	
	//Gaussian for y
	normal_distribution<double> dist_y(y, std_y);
	
	//Gaussian for theta
	normal_distribution<double> dist_theta(theta, std_theta);
	
	//Init all the particles
	for(int i=0; i<num_particles; i++)
	{
	    Particle p_ini;
	    
	    p_ini.id = i;
	    p_ini.x = dist_x(gen);
	    p_ini.y = dist_y(gen);
	    p_ini.theta = dist_theta(gen);
	    p_ini.weight = 1;
	    
	    particles.push_back(p_ini);
	    weights.push_back(1);
	}   //end for
	
	is_initialized = true;

}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/
	
	//Random Gaussian Noise (zero mean)
	default_random_engine gen;
	double std_x = std_pos[0];
	double std_y = std_pos[1];
	double std_theta = std_pos[2];
	
	//Gaussian for x
	//normal_distribution<double> dist_x(0, std_x);
	
	//Gaussian for y
	//normal_distribution<double> dist_y(0, std_y);
	
	//Gaussian for theta
	//normal_distribution<double> dist_theta(0, std_theta);
		
	
	//Prediction next time step
	for(int i=0; i<particles.size(); i++)
	{	    
	    //Bycicle Model
	    if(fabs(yaw_rate) < 0.0001)
	    {
	        particles[i].x += (velocity * delta_t) * cos(particles[i].theta);
	        particles[i].y += (velocity * delta_t) * sin(particles[i].theta);      
	    }
	    else
	    {	        	        
	        particles[i].x += (velocity/yaw_rate) * (sin(particles[i].theta + yaw_rate*delta_t) - sin(particles[i].theta));
	        particles[i].y += (velocity/yaw_rate)  *(cos(particles[i].theta) - cos(particles[i].theta + yaw_rate*delta_t));
	        particles[i].theta += yaw_rate*delta_t;
	        
	    } //end if
	    

	    //Gaussian for x
	    normal_distribution<double> dist_x(particles[i].x, std_x);
	
	    //Gaussian for y
	    normal_distribution<double> dist_y(particles[i].y, std_y);
	
	    //Gaussian for theta
	    normal_distribution<double> dist_theta(particles[i].theta, std_theta);
	    
	    particles[i].x = dist_x(gen);
	    particles[i].y = dist_y(gen);
	    particles[i].theta = dist_theta(gen);
	    
	    
	} //end for i particles
	

}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.
	
	std::vector<LandmarkObs> obs_r;
	
	for(int i = 0; i < observations.size(); i++)
	{
	    //Min distance is big, so then if we find a smaller one we update it.
	    double min_distance = numeric_limits<double>::max();
	    
	    //cout << " OBSERVATION OLD " <<   observations[i].id<< endl; 
		    
	    for(int j = 0; j < predicted.size(); j++)
	    {        
	        double distance = dist(observations[i].x, observations[i].y, predicted[j].x, predicted[j].y); 
	        
	        //if more than one, I just take the last
	        if(distance < min_distance)
	        {
	           //set min distance (big) with new distance (smaller)
	           min_distance = distance;
	           
	           observations[i].id =  predicted[j].id;
	           
	        }   //end if
	        
	    }   //end for j predicted
  
	}   //end for i observations

}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
		const std::vector<LandmarkObs> &observations, const Map &map_landmarks) {
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
	
	//Convert observations to map coordinates
	
	//Particles coordinates
	double x_p = 0;
	double y_p = 0;
	double theta_p = 0;
	
    //Observations coordinates
	double x_c = 0;
	double y_c = 0;
	//double theta_c = 0;
	
	//Map coordinates
	double x_m = 0;
	double y_m = 0;
	//double theta_m = 0;
	
	//vector of transformed observations
	std::vector<LandmarkObs> obs_t;
	
	//convert map landmarks to landmarks
	std::vector<LandmarkObs> map_pred;
	
	//filtered map landmarks (those that are within the sensor range)
	std::vector<LandmarkObs> map_pred_filt;
	
	double sig_x = std_landmark[0];
	double sig_y = std_landmark[1];
	
	//normalization term
	double gauss_norm = (1/(2 * M_PI * sig_x * sig_y));
	
	
	for(int m=0; m<map_landmarks.landmark_list.size(); m++)
	{
	    Map::single_landmark_s landmark_map = map_landmarks.landmark_list[m];
	    LandmarkObs landmark_pred;
	    
	    landmark_pred.id = landmark_map.id_i;
	    landmark_pred.x = landmark_map.x_f;
	    landmark_pred.y = landmark_map.y_f;
	    
	    //cout << landmark_pred.id << " : " << landmark_pred.x << " :  " << landmark_pred.y << endl;
	    
	    map_pred.push_back(landmark_pred);

	} //end for m map landmarks
	
	for(int i=0; i<particles.size(); i++)
	{
	    obs_t.clear();
	    
	    x_p = particles[i].x;
	    y_p = particles[i].y;
	    theta_p = particles[i].theta;
	    
	    //Convert observations to map coordinates
	    for(int j = 0; j<observations.size(); j++)
	    {
	        LandmarkObs obs;
	        
	        x_c = observations[j].x;
	        y_c = observations[j].y;
	        
	          
	        //convert to map coordinates
	        x_m = x_p + (cos(theta_p) * x_c) - (sin(theta_p) * y_c);
	        y_m = y_p + (sin(theta_p) * x_c) + (cos(theta_p) * y_c);
	        
	        obs.id = observations[j].id;
	        obs.x = x_m;
	        obs.y = y_m; 
	        
	        obs_t.push_back(obs);

	    }  //end for j observations
	    
	    //Just compare with map landmarks in the sensor range
	    map_pred_filt.clear();
	    for(int d = 0; d< map_pred.size(); d++)
	    {        
	        double dist_part_land = dist(map_pred[d].x, map_pred[d].y, particles[i].x, particles[i].y);	        
	        
	        //if the distance in greate, I delete that landmark from the vector
	        if(dist_part_land < sensor_range)
	        {
	            //cout << " RANGE " << dist_part_land << endl;
	            map_pred_filt.push_back(map_pred[d]); 	
	        }    //end if
	        
	     }   //end for d map_pred
	
	   
	    //Perform dataAssociation for that particle
	    this -> dataAssociation(map_pred_filt, obs_t);
	    
	    //cout << " MAP SIZE: " << map_pred_filt.size() << endl;
	    
	    //Now that I have the data association, I have to do the weight probabilities...
	    
	    //Go through the observations
	    for (int o = 0; o<obs_t.size(); o++)
	    {
	    
	        double x_o, y_o, mu_x, mu_y;
	        
	        x_o = obs_t[o].x;
	        y_o = obs_t[o].y;
	        
	        //this is the id of the associated prediction
	        int assoc_pred = obs_t[o].id;
	        
	        //Now I go through the predictions to get mu
	        for (int p = 0; p<map_pred_filt.size(); p++)
	        {
	            if(assoc_pred == map_pred_filt[p].id)
	            {
	                //cout << assoc_pred << " " << map_pred_filt[p].id << endl;
	                mu_x = map_pred_filt[p].x;
	                mu_y = map_pred_filt[p].y;
	            }
	        } //end for p map_pred_filt
	        
	        //exponent
	        double exponent = (pow((x_o - mu_x),2))/(2 * pow(sig_x,2)) + (pow((y_o - mu_y),2))/(2 * pow(sig_y,2));
	        
	        //cout << "EXP: " << exp(-exponent) << endl;
	        //cout << "NORM: " << gauss_norm << endl;

	        double prob_w = gauss_norm * exp(-exponent);
	        
	        particles[i].weight *= prob_w;
	        weights[i] = particles[i].weight;
	    
	    } //end of o observations
	    
	    //cout << " PARTICLE: " << particles[i].id << " WEIGHT: " << weights[i] << " " << particles[i].weight << endl;
	    
	} //end for i particles
	
}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
	
	//vector<Particle> new_particles;
	//vector<double> weights;
	//std::discrete_distribution<int> ind();
	// Vector for new particles
	
	std::random_device rd_wts;
	std::mt19937 generator_wts(rd_wts());


	// Creates a discrete distribution for weight.
	std::discrete_distribution<int> distribution_wts(weights.begin(), weights.end());
	std::vector<Particle> resampled_particles;

	// Resample
	for(int i=0;i<num_particles;i++){
		Particle particles_i = particles[distribution_wts(generator_wts)];
		resampled_particles.push_back(particles_i);
	}
	particles = resampled_particles;
  
	

}

Particle ParticleFilter::SetAssociations(Particle particle, std::vector<int> associations, std::vector<double> sense_x, std::vector<double> sense_y)
{
	//particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
	// associations: The landmark id that goes along with each listed association
	// sense_x: the associations x mapping already converted to world coordinates
	// sense_y: the associations y mapping already converted to world coordinates

	//Clear the previous associations
	particle.associations.clear();
	particle.sense_x.clear();
	particle.sense_y.clear();

	particle.associations= associations;
 	particle.sense_x = sense_x;
 	particle.sense_y = sense_y;

 	return particle;
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
