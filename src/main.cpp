#include <math.h>
#include <uWS/uWS.h>
#include <chrono>
#include <iostream>
#include <thread>
#include <vector>
#include "Eigen-3.3/Eigen/Core"
#include "Eigen-3.3/Eigen/QR"
#include "MPC.h"
#include "json.hpp"

// for convenience
using json = nlohmann::json;

// Checks if the SocketIO event has JSON data.
// If there is data the JSON object in string format will be returned,
// else the empty string "" will be returned.
string hasData(const string &s) {
  auto found_null = s.find("null");
  auto b1 = s.find_first_of("[");
  auto b2 = s.rfind("}]");
  if (found_null != string::npos) {
    return "";
  }
  else if (b1 != string::npos && b2 != string::npos) {
    return s.substr(b1, b2 - b1 + 2);
  }

  return "";
}

int main() {
  // MPC parameters.
  int N = 8;
  double dt = 0.11;
  double vr = 40.0;

  // Initialize MPC to given parameters.
  MPC mpc(N, dt, vr);

  uWS::Hub h;
  h.onMessage([&mpc](uWS::WebSocket<uWS::SERVER> ws, char *data, size_t length, uWS::OpCode opCode) {
    // Receive status message from simulator
    // and log it to standard output.
    string sdata = string(data).substr(0, length);
    std::cout << sdata << std::endl;

    // "42" at the start of the message means there's a websocket message event.
    // The 4 signifies a websocket message
    // The 2 signifies a websocket event
    if (!(sdata.size() > 2 && sdata[0] == '4' && sdata[1] == '2')) {
      return;
    }

    string s = hasData(sdata);
    if (s == "") {
      return;
    }

    auto j = json::parse(s);
    string event = j[0].get<string>();
    if (event != "telemetry") {
      return;
    }

    // Current car state.
    double px = j[1]["x"];
    double py = j[1]["y"];
    double psi = j[1]["psi"];
    double vt = j[1]["speed"];

    // Current reference route.
    vector<double> ptsx = j[1]["ptsx"];
    vector<double> ptsy = j[1]["ptsy"];
    // find local bindus
    margbindu reference = ontrack(px, py, psi, {ptsx, ptsy});

    // Compute sequence of actuations to make car approach reference route,
    // then take just one pair and discard the rest.
    auto actuations = mpc(vt, reference);
    auto a = actuations[2]; // Get the second actuation command pair
    auto d = actuations[3]; // to account for the system's latency.

    // Compute throttle and steering angle values.
    json msgJson;
    double tt = j[1]["throttle"];
    msgJson["throttle"] = tt + 0.1 * a * mpc.dt;
    msgJson["steering_angle"] = -d;

    // Display reference margbindu.
    msgJson["next_x"] = reference[0];
    msgJson["next_y"] = reference[1];

    // Display planned margbindu.
    margbindu planned = perform(mpc.dt, vt, actuations);
    msgJson["mpc_x"] = planned[0];
    msgJson["mpc_y"] = planned[1];

    // Finish building steering message
    // and log it to standard output.
    auto msg = "42[\"steer\"," + msgJson.dump() + "]";
    std::cout << msg << std::endl;

    // Simulate system latency and send message.
    auto latency = chrono::milliseconds((int) (mpc.dt * 1000));
    this_thread::sleep_for(latency);
    ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
  });

  // We don't need this since we're not using HTTP but if it's removed the
  // program
  // doesn't compile :-(
  h.onHttpRequest([](uWS::HttpResponse *res, uWS::HttpRequest req, char *data, size_t, size_t) {
    const std::string s = "<h1>Hello world!</h1>";
    if (req.getUrl().valueLength == 1) {
      res->end(s.data(), s.length());
    }
    else {
      // i guess this should be done more gracefully?
      res->end(nullptr, 0);
    }
  });

  h.onConnection([&h](uWS::WebSocket<uWS::SERVER> ws, uWS::HttpRequest req) {
    std::cout << "Connected!!!" << std::endl;
  });

  h.onDisconnection([&h](uWS::WebSocket<uWS::SERVER> ws, int code, char *message, size_t length) {
    ws.close();
    std::cout << "Disconnected" << std::endl;
  });

  int port = 4567;
  if (h.listen(port)) {
    std::cout << "Listening to port " << port << std::endl;
  }
  else {
    std::cerr << "Failed to listen to port" << std::endl;
    return -1;
  }

  h.run();
}
