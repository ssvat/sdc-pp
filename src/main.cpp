#include <fstream>
#include <math.h>
#include <uWS/uWS.h>
#include <chrono>
#include <iostream>
#include <thread>
#include <vector>
#include "Eigen-3.3/Eigen/Core"
#include "Eigen-3.3/Eigen/QR"
#include "json.hpp"
#include "spline.h"  //Use Spline instead of polynomial fits. Spline.h is created separately.

using namespace std;

// for convenience
using json = nlohmann::json;

// For converting back and forth between radians and degrees.
constexpr double pi() { return M_PI; }
double deg2rad(double x) { return x * pi() / 180; }
double rad2deg(double x) { return x * 180 / pi(); }

// From uWebSocekt template to check if the SocketIO event has JSON data.
// If there is data the JSON object in string format will be returned,
// else the empty string "" will be returned.
string hasData(string s) {
  auto found_null = s.find("null");
  auto b1 = s.find_first_of("[");
  auto b2 = s.find_first_of("}");
  if (found_null != string::npos) {
    return "";
  } else if (b1 != string::npos && b2 != string::npos) {
    return s.substr(b1, b2 - b1 + 2);
  }
  return "";
}

//To calculate Euclidean distance between two points
double distance(double x1, double y1, double x2, double y2)
{
  return sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1));
}

// To see the closest waypoints on a map of all waypoints around highways
int ClosestWaypoint(double x, double y, const vector<double> &maps_x, const vector<double> &maps_y)
{

  double closestLen = 100000; //large number
  int closestWaypoint = 0;

  for(int i = 0; i < maps_x.size(); i++)
    {
      double map_x = maps_x[i];
      double map_y = maps_y[i];
      double dist = distance(x,y,map_x,map_y);
      if(dist < closestLen)
	{
	  closestLen = dist;
	  closestWaypoint = i;
	}

    }

  return closestWaypoint;

}

// To see the waypoints behind or in front of our car
int NextWaypoint(double x, double y, double theta, const vector<double> &maps_x, const vector<double> &maps_y)
{

  int closestWaypoint = ClosestWaypoint(x,y,maps_x,maps_y);

  double map_x = maps_x[closestWaypoint];
  double map_y = maps_y[closestWaypoint];

  double heading = atan2( (map_y-y),(map_x-x) );

  double angle = abs(theta-heading);
  // Estimate car angle where the car is looking at.
  if(angle > pi()/4)
    {
      closestWaypoint++;
    }

  return closestWaypoint;

}

// Transform from Cartesian x,y coordinates to Frenet s,d coordinates
vector<double> getFrenet(double x, double y, double theta, const vector<double> &maps_x, const vector<double> &maps_y)
{
  int next_wp = NextWaypoint(x,y, theta, maps_x,maps_y);

  int prev_wp;
  prev_wp = next_wp-1;
  if(next_wp == 0)
    {
      prev_wp  = maps_x.size()-1;
    }

  double n_x = maps_x[next_wp]-maps_x[prev_wp];
  double n_y = maps_y[next_wp]-maps_y[prev_wp];
  double x_x = x - maps_x[prev_wp];
  double x_y = y - maps_y[prev_wp];

  // find the projection of x onto n.
  // Each waypoint has an (x,y) global map position, a Frenet s value and Frent d unit normal vector.
  double proj_norm = (x_x*n_x+x_y*n_y)/(n_x*n_x+n_y*n_y);
  double proj_x = proj_norm*n_x;
  double proj_y = proj_norm*n_y;

  double frenet_d = distance(x_x,x_y,proj_x,proj_y); // Including the x and y components 

  //see if d value is positive or negative by comparing it to a center point

  double center_x = 1000-maps_x[prev_wp];
  double center_y = 2000-maps_y[prev_wp];
  double centerToPos = distance(center_x,center_y,x_x,x_y);
  double centerToRef = distance(center_x,center_y,proj_x,proj_y);

  if(centerToPos <= centerToRef)
    {
      frenet_d *= -1; //To be used to calculate lane positions.
    }

  // To calculate s value. s refers to the distance along the direction of the road.
  double frenet_s = 0; // Starting point as zero
  for(int i = 0; i < prev_wp; i++)
    {
      frenet_s += distance(maps_x[i],maps_y[i],maps_x[i+1],maps_y[i+1]); // Include (i+1)th component for transitioning, otherwise sitting still
    }

  frenet_s += distance(0,0,proj_x,proj_y);

  return {frenet_s,frenet_d};

}

// Transform from Frenet s,d coordinates to Cartesian x,y
vector<double> getXY(double s, double d, const vector<double> &maps_s, const vector<double> &maps_x, const vector<double> &maps_y)
{
  int prev_wp = -1;

  while(s > maps_s[prev_wp+1] && (prev_wp < (int)(maps_s.size()-1) ))
    {
      prev_wp++;
    }

  int wp2 = (prev_wp+1)%maps_x.size();

  double heading = atan2((maps_y[wp2]-maps_y[prev_wp]),(maps_x[wp2]-maps_x[prev_wp]));
  // the x,y,s along the segment
  double seg_s = (s-maps_s[prev_wp]);

  double seg_x = maps_x[prev_wp]+seg_s*cos(heading);
  double seg_y = maps_y[prev_wp]+seg_s*sin(heading);

  double perp_heading = heading-pi()/2;

  double x = seg_x + d*cos(perp_heading);
  double y = seg_y + d*sin(perp_heading);

  return {x,y};

}

int main() {
  uWS::Hub h;

  // Load up map values for waypoint's x,y,s and d normalized normal vectors
  vector<double> map_waypoints_x;
  vector<double> map_waypoints_y;
  vector<double> map_waypoints_s;
  vector<double> map_waypoints_dx;
  vector<double> map_waypoints_dy;

  // Input a waypoint map
  string map_file_ = "../data/highway_map.csv";
  // The max s value before wrapping around the track back to 0
  double max_s = 6945.554;
  //Load the map
  ifstream in_map_(map_file_.c_str(), ifstream::in);

  string line;
  while (getline(in_map_, line)) {
    istringstream iss(line);
    double x;
    double y;
    float s;
    float d_x;  //Normal components of waypoints
    float d_y;
    iss >> x;
    iss >> y;
    iss >> s;
    iss >> d_x;
    iss >> d_y;
    map_waypoints_x.push_back(x);
    map_waypoints_y.push_back(y);
    map_waypoints_s.push_back(s);
    map_waypoints_dx.push_back(d_x);
    map_waypoints_dy.push_back(d_y);
  }

  // Car's lane. Stating at middle lane. 1 refers to central lane.
  int lane = 1;

  // Reference velocity.
  double ref_vel = 0.0; // mph, starting at zero to avoid cold start

  h.onMessage([&ref_vel, &lane, &map_waypoints_x,&map_waypoints_y,&map_waypoints_s,&map_waypoints_dx,&map_waypoints_dy]
	      (uWS::WebSocket<uWS::SERVER> ws, char *data, size_t length, uWS::OpCode opCode) {



		// "42" at the start of the message means there's a websocket message event.
		// The 4 signifies a websocket message
		// The 2 signifies a websocket event
		//auto sdata = string(data).substr(0, length);
		//cout << sdata << endl;
		if (length && length > 2 && data[0] == '4' && data[1] == '2') {

		  auto s = hasData(data);

		  if (s != "") {
		    auto j = json::parse(s);  //json library

		    string event = j[0].get<string>();

		    if (event == "telemetry") {
		      // j[1] is the data JSON object

		      // Main car's localization Data
		      double car_x = j[1]["x"];
		      double car_y = j[1]["y"];
		      double car_s = j[1]["s"];
		      double car_d = j[1]["d"]; // Check if staying in the constant lane
		      double car_yaw = j[1]["yaw"];  //Car's angle
		      double car_speed = j[1]["speed"];
 
		      // Previous path data given to the Planner
		      auto previous_path_x = j[1]["previous_path_x"];
		      auto previous_path_y = j[1]["previous_path_y"];
		      // Previous path's end s and d values
		      double end_path_s = j[1]["end_path_s"];
		      double end_path_d = j[1]["end_path_d"];

		      // Sensor Fusion Data, a list of all other cars on the same side of the road (right-hand side).
		      auto sensor_fusion = j[1]["sensor_fusion"];

		      // Project begin:
		      // Define a path made up of (x,y) points that the car will visit sequentially every .02 seconds.
		      
		      // Set up size of previous path vector to help transition
		      int prev_size = previous_path_x.size();

		      // Preventing from collisions.
		      if (prev_size > 0) {
			car_s = end_path_s;
		      }

		      // Prediction : Analysing other cars positions.
		      bool car_ahead = false;
		      bool car_left = false;
		      bool car_righ = false;
		      //Find ref_v to use.
		      for ( int i = 0; i < sensor_fusion.size(); i++ ) {
			float d = sensor_fusion[i][6]; //d in Frenet coordinates
			int car_lane = -1;
			// If Car is on our lane.
			if ( d > 0 && d < 4 ) {
			  car_lane = 0;
			} else if ( d > 4 && d < 8 ) {
			  car_lane = 1;
			} else if ( d > 8 && d < 12 ) {
			  car_lane = 2;
			}
			if (car_lane < 0) {
			  continue;
			}
			// Find car speed.
			double vx = sensor_fusion[i][3]; //3rd element of velocity.
			double vy = sensor_fusion[i][4]; //4th element of velocity.
			double check_speed = sqrt(vx*vx + vy*vy);
			double check_car_s = sensor_fusion[i][5]; // s in Frenet coordinates
			// Estimate/Project car's future position after executing previous trajectory.
			// If using previous poitns, s value can be projected out to predcit car's future path.
			// Predicted Frenet s value will be its current s value plus its tranformed velocity multiplied by time elapsed.
			check_car_s += ((double)prev_size*0.02*check_speed);  //0.02 refers to 50 path points per second.

			if ( car_lane == lane ) {
			  // Car ahead and in our lane. Check s values greater than mine and s gap.
			  car_ahead |= check_car_s > car_s && check_car_s - car_s < 30;
			} else if ( car_lane - lane == -1 ) {
			  // Car at left
			  car_left |= car_s - 30 < check_car_s && car_s + 30 > check_car_s;
			} else if ( car_lane - lane == 1 ) {
			  // Car at right
			  car_righ |= car_s - 30 < check_car_s && car_s + 30 > check_car_s;
			}
		      }

		      // Behavior planning
		      double speed_diff = 0;
		      const double MAX_SPEED = 49.5;
		      const double MAX_ACC = .224;  //5 m/s2
		      if ( car_ahead ) { // Car ahead
			if ( !car_left && lane > 0 ) {
			  // if there is no car left and there is a left lane.
			  lane--; // Change lane left.
			} else if ( !car_righ && lane != 2 ){
			  // if there is no car right and there is a right lane.
			  lane++; // Change lane right.
			} else { // If getting too close to frount car, slow down by subtracting velocity increment.
			  speed_diff -= MAX_ACC;
			}
		      } else {
			if ( lane != 1 ) { // if we are not on the center lane.
			  if ( ( lane == 0 && !car_righ ) || ( lane == 2 && !car_left ) ) {
			    lane = 1; // Back to center.
			  }
			}
			if ( ref_vel < MAX_SPEED ) { // Starting acceleration slowly to avoid cold start.
			  speed_diff += MAX_ACC;
			}
		      }
		      // Constant velocity comes from the same distance increment being used.
		      // Create a list of widely spaced (x,y) waypoints, evenly spaced at 30 m in order to
		      // be incoporated with a spline and fill it in with more points that control speed.
		      vector<double> ptsx;
		      vector<double> ptsy;
                      // Reference x, y, and yaw states
		      // Either reference the starting points as where the car is or at the previous paths end point.
		      double ref_x = car_x;
		      double ref_y = car_y;
		      double ref_yaw = deg2rad(car_yaw);
                      // If previous size is almost empty, use the car as starting reference.  
		      if ( prev_size < 2 ) {
			// Calculate path points by using two points that make the path trangent to angle of car.
			double prev_car_x = car_x - cos(car_yaw);
			double prev_car_y = car_y - sin(car_yaw);

			ptsx.push_back(prev_car_x);
			ptsx.push_back(car_x);

			ptsy.push_back(prev_car_y);
			ptsy.push_back(car_y);
		      } else {
			// Use the last two points. Redefine reference state as previous path end point.
			ref_x = previous_path_x[prev_size - 1]; // Last point
			ref_y = previous_path_y[prev_size - 1];

			double ref_x_prev = previous_path_x[prev_size - 2]; // Last second point
			double ref_y_prev = previous_path_y[prev_size - 2];
			ref_yaw = atan2(ref_y-ref_y_prev, ref_x-ref_x_prev);
                        // Use two pionts that make the path tangent to the previous path's end point
			ptsx.push_back(ref_x_prev);
			ptsx.push_back(ref_x);

			ptsy.push_back(ref_y_prev);
			ptsy.push_back(ref_y);
		      }

		      // In Frenet coordinates, add evenly 30 m spaced points ahead of the starting reference to stay on the track.
		      // Adding spaces makes car lane shift smooth. Also use getXY to stay in a lane.
		      // Lane width is 4 so another lane center would be 2 + 4*lane.
		      // In total, there are five points (3 here plus previous 2 reference points).
		      vector<double> next_wp0 = getXY(car_s + 30, 2 + 4*lane, map_waypoints_s, map_waypoints_x, map_waypoints_y);
		      vector<double> next_wp1 = getXY(car_s + 60, 2 + 4*lane, map_waypoints_s, map_waypoints_x, map_waypoints_y);
		      vector<double> next_wp2 = getXY(car_s + 90, 2 + 4*lane, map_waypoints_s, map_waypoints_x, map_waypoints_y);
                      // Boundary points
		      ptsx.push_back(next_wp0[0]); // First element
		      ptsx.push_back(next_wp1[0]);
		      ptsx.push_back(next_wp2[0]);

		      ptsy.push_back(next_wp0[1]); // Second element
		      ptsy.push_back(next_wp1[1]);
		      ptsy.push_back(next_wp2[1]);

		      // Shift the car's reference angle to local car coordinates, resulting a rotation angle of 0.
		      // Shift first and then rotation
		      for ( int i = 0; i < ptsx.size(); i++ ) {
			double shift_x = ptsx[i] - ref_x;
			double shift_y = ptsy[i] - ref_y;

			ptsx[i] = shift_x * cos(0 - ref_yaw) - shift_y * sin(0 - ref_yaw);
			ptsy[i] = shift_x * sin(0 - ref_yaw) + shift_y * cos(0 - ref_yaw);
		      }

		      // Create the spline to interpolate points.
		      tk::spline s;
		      // Set points x and points y to the spline. Add those 5 anchor points.
		      s.set_points(ptsx, ptsy);

                      // To get the car going by defining the actual (x,y) points used in the planner.
		      vector<double> next_x_vals;
		      vector<double> next_y_vals;
		      
		      // Add a new element at the end of the path vector, after last element of previous path.
		      // Loop and start with previous pass poitns from last time, preventing from a spike in accleration.
		      for ( int i = 0; i < prev_size; i++ ) {
			next_x_vals.push_back(previous_path_x[i]);  // To go in a straight line continuously
			next_y_vals.push_back(previous_path_y[i]);
		      }

		      // Calculate how to break up spline points in order to travle at the desired reference velocity.
		      double target_x = 30.0; //Horizontal distance
		      double target_y = s(target_x); // Ask the spline what's the y for the given x
		      //Distance between car itslef to that target.
		      double target_dist = sqrt(target_x*target_x + target_y*target_y); 

		      double x_add_on = 0; //For local transformation.
                      // Only loop for current path points, excluding simulated, previous path points.
		      // There are 50 path points in total. Fill up the rest of the path planner after filling it with previous points

		      for( int i = 1; i < 50 - prev_size; i++ ) {
			ref_vel += speed_diff;
			if ( ref_vel > MAX_SPEED ) {
			  ref_vel = MAX_SPEED;
			} else if ( ref_vel < MAX_ACC ) {
			  ref_vel = MAX_ACC;
			}
			// More efficient by considering velocity increment
			double N = target_dist/(0.02*ref_vel/2.24); 
			double x_point = x_add_on + target_x/N;
			double y_point = s(x_point);

			x_add_on = x_point;

			double x_ref = x_point;
			double y_ref = y_point;
			
                        // Rotation back to global coordinates
			x_point = x_ref * cos(ref_yaw) - y_ref * sin(ref_yaw);
			y_point = x_ref * sin(ref_yaw) + y_ref * cos(ref_yaw);

			x_point += ref_x;
			y_point += ref_y;
			// Append new waypoints, until the new path has 50 total waypoints.
			next_x_vals.push_back(x_point);
			next_y_vals.push_back(y_point);
		      }

		      json msgJson;  //Library
		      
                      // Define a path made up of (x,y) points that the car will visit sequentially every .02 seconds.
		      
		      msgJson["next_x"] = next_x_vals;
		      msgJson["next_y"] = next_y_vals;

		      auto msg = "42[\"control\","+ msgJson.dump()+"]";

		      //this_thread::sleep_for(chrono::milliseconds(1000));
		      ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);

		    }
		  } else {
		    // Manual driving
		    std::string msg = "42[\"manual\",{}]";
		    ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
		  }
		}
	      });

  // We don't need this since we're not using HTTP but if it's removed the
  // program
  // doesn't compile :-(
  h.onHttpRequest([](uWS::HttpResponse *res, uWS::HttpRequest req, char *data,
		     size_t, size_t) {
		    const std::string s = "<h1>Hello world!</h1>";
		    if (req.getUrl().valueLength == 1) {
		      res->end(s.data(), s.length());
		    } else {
		      // i guess this should be done more gracefully?
		      res->end(nullptr, 0);
		    }
		  });

  h.onConnection([&h](uWS::WebSocket<uWS::SERVER> ws, uWS::HttpRequest req) {
      std::cout << "Connected!!!" << std::endl;
    });

  h.onDisconnection([&h](uWS::WebSocket<uWS::SERVER> ws, int code,
			 char *message, size_t length) {
		      ws.close();
		      std::cout << "Disconnected" << std::endl;
		    });

  int port = 4567;
  if (h.listen(port)) {
    std::cout << "Listening to port " << port << std::endl;
  } else {
    std::cerr << "Failed to listen to port" << std::endl;
    return -1;
  }
  h.run();
}
